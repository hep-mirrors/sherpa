#include "AddOns/OpenLoops/EWVirtKFactor_Setter.H"

#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Flavour.H"

#include "MODEL/Main/Model_Base.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"


using namespace OpenLoops;
using namespace PHASIC;
using namespace ATOOLS;

EWVirtKFactor_Setter::EWVirtKFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args), p_ewloop(NULL), m_deltaew(0.)
{
  DEBUG_FUNC("");
  InitEWVirt();
}

EWVirtKFactor_Setter::~EWVirtKFactor_Setter()
{
  if (p_ewloop) { delete p_ewloop; p_ewloop=NULL; }
}

double EWVirtKFactor_Setter::KFactor()
{
  if(!m_on) return 1.;
  CopyMomenta();
  CalcEWCorrection();
  return (1.+m_deltaew);
}

double EWVirtKFactor_Setter::KFactor(const NLO_subevt& evt)
{
  if(!m_on) return 1.;
  CopyMomenta(evt);
  CalcEWCorrection();
  return (1.+m_deltaew);
}

void EWVirtKFactor_Setter::CopyMomenta()
{
  m_p=p_proc->Integrator()->Momenta();
}

void EWVirtKFactor_Setter::CopyMomenta(const NLO_subevt& evt)
{
  Vec4D_Vector p(evt.p_mom, &(evt.p_mom[evt.m_n]));
  for (size_t i(0);i<p_proc->NIn();++i) p[i]=-p[i];
}

void EWVirtKFactor_Setter::InitEWVirt()
{
  Process_Info loop_pi(p_proc->Info());
  if (p_proc->Info().m_fi.m_nloewtype!=nlo_type::lo ||
      !(p_proc->Info().m_fi.m_nloqcdtype==nlo_type::lo ||
        p_proc->Info().m_fi.m_nloqcdtype&nlo_type::born)) return;
  loop_pi.m_fi.m_nloewtype=nlo_type::loop;
  loop_pi.m_loopgenerator="OpenLoops";
  loop_pi.m_maxcpl=p_proc->MaxOrders();
  loop_pi.m_mincpl=p_proc->MinOrders();
  ++loop_pi.m_maxcpl[1];
  ++loop_pi.m_mincpl[1];
  msg_Debugging()<<"Load OpenLoops process for "<<p_proc->Name()
                 <<" of order "<<loop_pi.m_mincpl<<" .. "<<loop_pi.m_maxcpl
                 <<std::endl;
  p_ewloop=PHASIC::Virtual_ME2_Base::GetME2(loop_pi);
  if (!p_ewloop) {
    THROW(not_implemented,"Couldn't find OpenLoops EW Virtual for "
                          +p_proc->Name());
  }
  MODEL::s_model->GetCouplings(m_cpls);
  p_ewloop->SetCouplings(m_cpls);
}

void EWVirtKFactor_Setter::CalcEWCorrection()
{
  DEBUG_FUNC("");
  m_deltaew=0.;
  p_ewloop->SetRenScale(p_proc->ScaleSetter()->Scale(stp::ren,1));
  p_ewloop->Calc(m_p);
  // OL returns V/(as/2pi*B)
  double fac(p_ewloop->AlphaQCD()/2./M_PI);
  double B(p_ewloop->ME_Born());
  double V(p_ewloop->ME_Finite()*fac*B);
  m_deltaew=V/B;
  if (msg_LevelIsDebugging()) {
    msg_Out()<<" p_T    = "<<(m_p[2]+m_p[3]).PPerp()<<std::endl;
    msg_Out()<<" \\mu_R  = "<<p_proc->ScaleSetter()->Scale(stp::ren,1)<<std::endl;
    msg_Out()<<" VI_e2  = "<<p_ewloop->ME_E2()*fac*B<<std::endl;
    msg_Out()<<" VI_e1  = "<<p_ewloop->ME_E1()*fac*B<<std::endl;
    msg_Out()<<" VI_fin = "<<V<<std::endl;
    msg_Out()<<" B      = "<<B<<std::endl;
    msg_Out()<<" \\delta = "<<m_deltaew<<std::endl;
  }
}


DECLARE_GETTER(EWVirtKFactor_Setter,"EWVirt",
               KFactor_Setter_Base,
               KFactor_Setter_Arguments);

KFactor_Setter_Base *Getter<KFactor_Setter_Base,
                            KFactor_Setter_Arguments,
                            EWVirtKFactor_Setter>::operator()
(const KFactor_Setter_Arguments &args) const
{
  msg_Info()<<"Loading EWVirt KFactor for "<<args.p_proc->Name()<<std::endl;
  return new EWVirtKFactor_Setter(args);
}

void Getter<KFactor_Setter_Base,
            KFactor_Setter_Arguments,
            EWVirtKFactor_Setter>::PrintInfo(std::ostream &str,
                                             const size_t width) const
{
  str<<"EWVirt K-Factor\n";
}
