#include "AddOns/GoSam/GoSam_Virtual.H"

#include "AddOns/GoSam/GoSam_Interface.H"
#include "ATOOLS/Math/Poincare.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace GoSam;
using namespace std;


GoSam_Virtual::GoSam_Virtual(const Process_Info& pi,
                             const Flavour_Vector& flavs,bool active,
                             int proc_id) :
  Virtual_ME2_Base(pi, flavs), m_pn(flavs.size()), m_active(active),
  m_proc_id(proc_id)
{
  m_Norm = pi.m_ii.ISSymmetryFactor()*pi.m_fi.FSSymmetryFactor();
  m_oqcd = pi.m_maxcpl[0];
  m_oqed = pi.m_maxcpl[1];
  if (!m_active) return;
  m_mode = 0;


  m_g2qed = 4.0*M_PI*MODEL::s_model->ScalarConstant("alpha_QED");
  m_g2qcd = 4.0*M_PI*MODEL::s_model->ScalarConstant("alpha_S");

  m_momenta = vector<double>(m_pn*5);
  for (size_t i=0;i<m_pn;i++) m_momenta[4+i*5]=flavs[i].Mass();
  m_result = {0.,0.,0.,1.};

  auto s = Settings::GetMainSettings()["GOSAM"];
  string irr = s["IRregularisation"].Get<string>();
  if(irr == "CDR") m_drmode = 0;
  else if(irr == "DRED") m_drmode = 1;

#ifdef VALIDATE_GOSAM_INTERFACE
  m_errorOnFailedCheck = s["ErrorOnFailedCheck"].SetDefault(false).Get<bool>();
  m_printAllChecks = s["PrintAllChecks"].SetDefault(false).Get<bool>();
#endif

}

void GoSam_Virtual::Calc(const Vec4D_Vector& pp) {
  DEBUG_FUNC(m_proc_id);
  if (!m_active) return;
  if (m_proc_id<0) return;

  // TODO: Do we need CMS boost?
  Vec4D_Vector momenta(pp);
  Poincare cms(momenta[0]+momenta[1]);
  for (size_t i(0);i<momenta.size();++i) {
    msg_Debugging()<<momenta[i]<<"\n";
    cms.Boost(momenta[i]);
  }

  for (size_t i=0;i<m_pn;i++) {
    m_momenta[0+i*5]=momenta[i][0];
    m_momenta[1+i*5]=momenta[i][1];
    m_momenta[2+i*5]=momenta[i][2];
    m_momenta[3+i*5]=momenta[i][3];
  }
  m_result = {0.,0.,0.,0.};
  double scale = sqrt(m_mur2);
  double acc = 0;
  double ierr = 1;

  GoSam::GoSam_Interface::SetParameter("alphaS", AlphaQCD(), ierr);
  GoSam::GoSam_Interface::SetParameter("alpha", AlphaQED(), ierr);
  GoSam::OLP_EvalSubProcess2(&m_proc_id, m_momenta.data(), &scale, m_result.data(), &acc);
  m_born         = m_result[3]/(2.*M_PI)*AlphaQCD();
  m_res.Finite() = m_result[2]/m_born;
  m_res.IR()     = m_result[1]/m_born;
  m_res.IR2()    = m_result[0]/m_born;
  m_born         *= 2.*M_PI/AlphaQCD();

#ifdef VALIDATE_GOSAM_INTERFACE
  if(p_check) {
    MODEL::Coupling_Map cpls;
    cpls.insert({"Alpha_QCD",p_aqcd});
    cpls.insert({"Alpha_QED",p_aqed});
    p_check->SetCouplings(cpls);
    p_check->SetRenScale(m_mur2);
    p_check->Calc(momenta);
    bool checkfailed=false;
    if(!IsEqual(m_res.IR2(),p_check->ME_E2(),1e-2)) {
      msg_Out()<<"Check double pole GoSam ("<<m_res.IR2()
               <<") vs. Internal ("<<p_check->ME_E2()<<") failed.\n";
      checkfailed=true;
    }
    else if(m_printAllChecks) {
      msg_Out()<<"Check double pole GoSam ("<<m_res.IR2()<<") vs. Internal ("
               <<p_check->ME_E2()<<") succesful.\n";
    }

    if(!IsEqual(m_res.IR(),p_check->ME_E1(),1e-2)) {
      msg_Out()<<"Check single pole GoSam ("<<m_res.IR()
               <<") vs. Internal ("<<p_check->ME_E1()<<") failed.\n";
      checkfailed=true;
    }
    else if(m_printAllChecks) {
      msg_Out()<<"Check single pole GoSam ("<<m_res.IR()<<") vs. Internal ("
               <<p_check->ME_E1()<<") succesful.\n";
    }

    if(!IsEqual(m_res.Finite(),p_check->ME_Finite(),1e-2)) {
      msg_Out()<<"Check finite loop GoSam ("<<m_res.Finite()
               <<") vs. Internal ("<<p_check->ME_Finite()<<") failed.\n";
      checkfailed=true;
    }
    else if(m_printAllChecks) {
      msg_Out()<<"Check finite loop GoSam ("<<m_res.Finite()<<") vs. Internal ("
               <<p_check->ME_Finite()<<") succesful.\n";
    }


    if(!IsEqual(m_born,p_check->ME_Born(),1e-2)) {
      msg_Out()<<"Check born GoSam ("<<m_born<<") vs. Internal ("
               <<p_check->ME_Born()<<") failed.\n";
      // No Born in (some) internal MEs.
      // checkfailed=true;
    }
    else if(m_printAllChecks) {
      msg_Out()<<"Check born GoSam ("<<m_born<<") vs. Internal ("
               <<p_check->ME_Born()<<") succesful.\n";
    }

    if(m_errorOnFailedCheck && checkfailed)
      THROW(fatal_error,"Failed check is not acceptable.");
  }
#endif

}

bool GoSam_Virtual::SetColours(const ATOOLS::Vec4D_Vector& momenta) {
  return true;
}

double GoSam_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
  return 4.*M_PI;
}

DECLARE_VIRTUALME2_GETTER(GoSam_Virtual,"GoSam")
Virtual_ME2_Base *ATOOLS::Getter<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,
				 GoSam_Virtual>::
operator()(const PHASIC::Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="GoSam") return nullptr;
  if (!(pi.m_fi.m_nlotype&nlo_type::loop)) return nullptr;
  Flavour_Vector fl=pi.ExtractFlavours();
  GoSam_Virtual * ret = new GoSam_Virtual(pi, fl, true, GoSam_Interface::RegisterProcess(pi,fl,true));
#ifdef VALIDATE_GOSAM_INTERFACE
  PHASIC::Process_Info picheck(pi);
  picheck.m_loopgenerator="Internal";
  ret->p_check = Virtual_ME2_Base::GetME2(picheck);
  if(ret->p_check)
    msg_Out()<<"Initialized Internal ME to validate GoSam interface.\n";
  else
    msg_Out()<<"Failed to initialize Internal ME to validate GoSam interface.\n";
#endif
  if(ret) {
    msg_Debugging()<<"Return GoSam Virtual.\n";
  }
  return ret;
}
