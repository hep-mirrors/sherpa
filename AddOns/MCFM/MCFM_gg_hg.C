#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  // README:
  // For Higgs production, choose model: MODEL = SM+EHC
  // It is important for the Higgs production to have all five flavours 
  // in the initial state, but the Yukawa coupling of the b must be
  // switched off:  YUKAWA[5] = 0.
  // Also, MCFM acts in the limit of mt->infinity,
  // thus a further correction term has been introduced
  //
  // We could actually also extend this to BSM models.


  class MCFM_gg_hg: public PHASIC::Virtual_ME2_Base {
  private:
    int                     m_pID;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    double                  m_mh2,m_Gh2,m_mZ2,m_GZ2,m_mW2,m_GW2;
    double                  m_ehcscale2,m_cplcorr,m_normcorr;


    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_gg_hg(const int & pID,
	      const PHASIC::Process_Info& pi,
	      const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_gg_hg();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void gg_hg_v_(double *p,double *msqv); 
  void gg_hwwg_v_(double *p,double *msqv); 
  void gg_hzzg_v_(double *p,double *msqv); 
  void gg_hgagag_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_gg_hg::MCFM_gg_hg(const int & pID,const Process_Info& pi,
		     const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  m_mh2(sqr(Flavour(kf_h0).Mass())),
  m_Gh2(sqr(Flavour(kf_h0).Width())),
  m_mZ2(sqr(Flavour(kf_Z).Mass())),
  m_GZ2(sqr(Flavour(kf_Z).Width())),
  m_mW2(sqr(Flavour(kf_Wplus).Mass())),
  m_GW2(sqr(Flavour(kf_Wplus).Width())),
  m_ehcscale2(MODEL::s_model->ScalarConstant(std::string("EHC_SCALE2"))),
  m_cplcorr(ewcouple_.vevsq/
	    sqr(MODEL::s_model->ScalarConstant(std::string("vev")))*
            sqr((*p_as)(m_ehcscale2)/qcdcouple_.as)*
	    sqr(MODEL::s_model->ScalarConstant(std::string("h0_gg_fac"))/
		(2./3.))),
  m_normcorr(4.*9./qcdcouple_.ason2pi)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  switch (m_pID) {
  case 204:
    m_cplcorr *= 
      sqr(Flavour(kf_tau).Yuk())/masses_.mbsq/3. *
      4.*sqr(masses_.wmass)/ewcouple_.gwsq/
      sqr(MODEL::s_model->ScalarConstant(std::string("vev")));
    break;
  case 208:
  case 209:
    m_cplcorr *=
      pow(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_QED"))/
	  MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))/
	  ewcouple_.gwsq,3.);
    break;
  case 210:
    m_cplcorr *= 
      sqr(MODEL::s_model->ScalarFunction(std::string("alpha_QED"),m_ehcscale2)/
	  (ewcouple_.esq/(4.*M_PI))) *
      1./(sqrt(2.)*ewcouple_.Gf*
	  sqr(MODEL::s_model->ScalarConstant(std::string("vev"))));
    break;
  }

  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_gg_hg::~MCFM_gg_hg()
{
  delete [] p_p;
  delete [] p_msqv;
}


double MCFM_gg_hg::CallMCFM(const int & i,const int & j) {
  //msg_Out()<<METHOD<<"("<<i<<", "<<j<<") for "<<m_pID<<".\n";
  switch (m_pID) {
  case 204: gg_hg_v_(p_p,p_msqv); break;
  case 208: gg_hwwg_v_(p_p,p_msqv); break;
  case 209: gg_hzzg_v_(p_p,p_msqv); break;
  case 210: gg_hgagag_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_gg_hg::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr*(*p_as)(m_mur2)/qcdcouple_.as);
  double sh((m_pID==204)?
	    (p[2]+p[3]).Abs2() :
	    (p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= 
    (sqr(sh-sqr(masses_.hmass))+
     sqr(masses_.hmass*masses_.hwidth))/
    (sqr(sh-m_mh2)+m_mh2*m_Gh2);
  if (m_pID==208 || m_pID==209) {
    double s1((p[2]+p[3]).Abs2()),s2((p[4]+p[5]).Abs2());
    if (m_pID==208) {
      corrfactor *= 
	(sqr(s1-sqr(masses_.wmass))+
	 sqr(masses_.wmass*masses_.wwidth))/
	(sqr(s1-m_mW2)+m_mW2*m_GW2);
      corrfactor *= 
	(sqr(s2-sqr(masses_.wmass))+
	 sqr(masses_.wmass*masses_.wwidth))/
	(sqr(s2-m_mW2)+m_mW2*m_GW2);
    }
    else {
      corrfactor *= 
	(sqr(s1-sqr(masses_.zmass))+
	 sqr(masses_.zmass*masses_.zwidth))/
	(sqr(s1-m_mZ2)+m_mZ2*m_GZ2);
      corrfactor *= 
	(sqr(s2-sqr(masses_.zmass))+
	 sqr(masses_.zmass*masses_.zwidth))/
	(sqr(s2-m_mZ2)+m_mZ2*m_GZ2);
    }
  }
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; corrfactor *= 8./3.; }
  if (j==21) { j=0; corrfactor *= 8./3.; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  epinv_.epinv=epinv2_.epinv2=0.0;
  double res(CallMCFM(i,j)  * corrfactor);
  epinv_.epinv=1.0;
  double res1(CallMCFM(i,j) * corrfactor);
  epinv2_.epinv2=1.0;
  double res2(CallMCFM(i,j) * corrfactor);
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);

  //msg_Out()<<METHOD<<" yields "<<m_res.Finite()
  //	   <<" + 1/eps * "<<m_res.IR()
  //	   <<" + 1/eps^2 * "<<m_res.IR2()
  //	   <<" for mb = "<<masses_.mb
  //	   <<" and "<<nflav_.nflav<<" active flavours"
  //	   <<" in "<<scheme_.scheme<<" ...  .\n";
}

double MCFM_gg_hg::Eps_Scheme_Factor(const Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gg_hg_Getter,"MCFM_gg_hg")
Virtual_ME2_Base *MCFM_gg_hg_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  else msg_Out()<<".\n";
  if (MODEL::s_model->Name()!=std::string("SM+EHC") ||
      Flavour(kf_b).Yuk()>0. ||
      !Flavour(kf_h0).IsOn())                           return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (!(pi.m_fi.m_nloqcdtype&nlo_type::loop))           return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (!(fl[0].Strong() && fl[1].Strong()))              return NULL;
  if (pi.m_fi.m_ps.size()!=2)                           return NULL;
  Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
  if (!flh==Flavour(kf_h0))                             return NULL;

  msg_Out()<<"Check numbers: "
	   <<fl.size()<<" external particles, "
	   <<pi.m_fi.m_ps.size()<<" props";
  if (pi.m_fi.m_ps.size()==0) msg_Out()<<".\n";
  else {
    msg_Out()<<":\n";
    for (size_t i=0;i<pi.m_fi.m_ps.size();i++) {
      msg_Out()<<"     "<<pi.m_fi.m_ps[i].m_ps.size()
	       <<" final state particles for propagator "<<i<<" "
	       <<"with flavour "<<pi.m_fi.m_ps[i].m_fl[0]<<".\n";
      for (size_t j=0;j<pi.m_fi.m_ps[i].m_ps.size();j++) {
	msg_Out()<<"     "<<pi.m_fi.m_ps[i].m_ps[j].m_ps.size()
		 <<" final state particles for propagator "<<i<<"["<<j<<"] "
		 <<"with flavour "<<pi.m_fi.m_ps[i].m_ps[j].m_fl[0]<<".\n";
      }
    }
  }

  int pID(0);
  //////////////////////////////////////////////////////////////////
  // tau tau final states
  //////////////////////////////////////////////////////////////////
  if (fl.size()==5 && pi.m_fi.m_ps.size()==2 && fl[4].Strong() &&
      fl[2]==fl[3].Bar() && fl[2].Kfcode()==15) {
    if (Flavour(kf_tau).Yuk()<=0.) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Setup for gg->[h->tau tau] (+jet), but tau Yukawa = 0."
		 <<std::endl;
      THROW(fatal_error,"Inconsistent setup.");
    }
    pID = 204;
  }
  //////////////////////////////////////////////////////////////////
  // photon photon final states 
  //////////////////////////////////////////////////////////////////
  else if (fl.size()==5 && pi.m_fi.m_ps.size()==2 && fl[4].Strong() &&
	   fl[2]==fl[3] && fl[2].Kfcode()==22) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   gg->[h->gamma gamma] not yet debugged.\n";
    THROW(fatal_error,"Not yet working.");
    if (fl.size()==4 && pi.m_fi.m_ps.size()==1)           pID = 210;
  }
  //////////////////////////////////////////////////////////////////
  // Z Z final states 
  //////////////////////////////////////////////////////////////////
  else if (fl.size()==7 && pi.m_fi.m_ps.size()==2 && fl[6].Strong() &&
	   pi.m_fi.m_ps[0].m_ps[0].m_fl[0].Kfcode()==23 &&
	   pi.m_fi.m_ps[0].m_ps[0].m_ps.size()==2 &&
	   pi.m_fi.m_ps[0].m_ps[1].m_fl[0].Kfcode()==23 &&
	   pi.m_fi.m_ps[0].m_ps[1].m_ps.size()==2) {
    Flavour Z11,Z12,Z21,Z22;
    Z11 = pi.m_fi.m_ps[0].m_ps[0].m_ps[0].m_fl[0];
    Z12 = pi.m_fi.m_ps[0].m_ps[0].m_ps[1].m_fl[0];
    Z21 = pi.m_fi.m_ps[0].m_ps[1].m_ps[0].m_fl[0];
    Z22 = pi.m_fi.m_ps[0].m_ps[1].m_ps[1].m_fl[0];
    msg_Out()<<"   test for ZZ: \n"
	     <<"   {"<<Z11<<", "<<Z12<<"} {"<<Z21<<", "<<Z22<<"}.\n";
    if (!(Z11==Z12.Bar() && Z21==Z22.Bar())) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   wrong flavour pairs: "
		 <<Z11<<" & "<<Z12<<", "<<Z21<<" & "<<Z22<<".\n"; 
      THROW(fatal_error,"wrong process.");
      return NULL;
    }
    if (Z11.Kfcode()==15 || Z21.Kfcode()==15) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   process should not involve taus: "
		 <<Z11<<" & "<<Z12<<", "<<Z21<<" & "<<Z22<<".\n"; 
      THROW(fatal_error,"wrong process.");
      return NULL;
    }
    if (Z11.Charge()!=0 && Z21.Charge()!=0)      pID = 209;
    else if (Z11.Charge()!=0 && Z21.Charge()==0) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   cannot deal with Z->nu nu yet.\n"
		 <<"   Must rescale couplings, will be done soon.\n";
      THROW(fatal_error,"wrong process.");
    }
  }  
  //////////////////////////////////////////////////////////////////
  // W W final states 
  //////////////////////////////////////////////////////////////////
  else if (fl.size()==7 && pi.m_fi.m_ps.size()==2 && fl[6].Strong() &&
	   pi.m_fi.m_ps[0].m_ps[0].m_fl[0].Kfcode()==24 &&
	   pi.m_fi.m_ps[0].m_ps[0].m_ps.size()==2 &&
	   pi.m_fi.m_ps[0].m_ps[1].m_fl[0].Kfcode()==24 &&
	   pi.m_fi.m_ps[0].m_ps[1].m_ps.size()==2) {
    Flavour W11,W12,W21,W22;
    W11 = pi.m_fi.m_ps[0].m_ps[0].m_ps[0].m_fl[0];
    W12 = pi.m_fi.m_ps[0].m_ps[0].m_ps[1].m_fl[0];
    W21 = pi.m_fi.m_ps[0].m_ps[1].m_ps[0].m_fl[0];
    W22 = pi.m_fi.m_ps[0].m_ps[1].m_ps[1].m_fl[0];
    msg_Out()<<"   test for WW: \n"
	     <<"   {"<<W11<<", "<<W12<<"} {"<<W21<<", "<<W22<<"}.\n";
    if (!W11.IsAnti() && W12.IsAnti() && !W21.IsAnti() && W22.IsAnti() &&
	((W11.Kfcode()==11 && W12.Kfcode()==12) || 
	 (W11.Kfcode()==13 && W12.Kfcode()==14)) && 
	((W21.Kfcode()==12 && W22.Kfcode()==11) || 
	 (W21.Kfcode()==14 && W22.Kfcode()==13)))        pID = 208;
  }
  //////////////////////////////////////////////////////////////////
  // check for 5 light flavours, this is hard-coded in MCFM
  //////////////////////////////////////////////////////////////////
  Flavour lq(kf_quark);
  double nlf(0.);
  for (size_t i(0); i<lq.Size(); ++i) if (!lq[i].IsMassive()) nlf++;
  nlf/=2.;
  if (nlf!=5.) {
    msg_Error()<<METHOD<<"(): MCFM only available for five light quark "
                       <<"flavours for this process"<<std::endl;
    THROW(fatal_error, "Five-flavour calculations with MCFM only.")
  }
  if (pID>0) {
    zerowidth_.zerowidth=true;
    if (nproc_.nproc>=0) {
      if (nproc_.nproc!=pID)
	THROW(not_implemented,
	      "Only one process class allowed when using MCFM");
    }
    nproc_.nproc=pID;
    chooser_();
    msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
    return new MCFM_gg_hg(pID,pi,fl);
  }
  return NULL;
}
