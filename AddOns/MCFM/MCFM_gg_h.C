#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  // README:
  // For Higgs production, choose model: MODEL = SM+EHC
  // It is important for the Higgs production to have all five flavours 
  // in the initial state, but the Yukawa coupling of the b must be
  // switched off:  YUKAWA_B = 0.  
  // Also, MCFM acts in the limit of mt->infinity,
  // thus a further correction term has been introduced
  //
  // We could actually also extend this to BSM models.


  class MCFM_gg_h: public PHASIC::Virtual_ME2_Base {
  private:
    int                     m_pID;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    double                  m_mh2,m_Gh2,m_cplcorr,m_normcorr;


    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_gg_h(const int & pID,
	      const PHASIC::Process_Info& pi,
	      const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_gg_h();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void gg_h_v_(double *p,double *msqv); 
  void qqb_hww_v_(double *p,double *msqv); 
  void qqb_hzz_v_(double *p,double *msqv); 
  void gg_hg_v_(double *p,double *msqv); 
  void gg_hwwg_v_(double *p,double *msqv); 
  void gg_hzzg_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_gg_h::MCFM_gg_h(const int & pID,const Process_Info& pi,
		     const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_cplcorr(ewcouple_.vevsq/
	    ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev")))*
	    ATOOLS::sqr((*p_as)(m_mh2)/qcdcouple_.as)),
  m_normcorr(4.*9./qcdcouple_.ason2pi)
{
  switch (m_pID) {
  case 112:
  case 204:
    m_cplcorr *= 
      ATOOLS::sqr(ATOOLS::Flavour(kf_tau).Yuk())/masses_.mbsq/3. *
      4.*ATOOLS::sqr(masses_.wmass)/ewcouple_.gwsq/
      ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev")));
    break;
  case 113:
  case 114:
  case 115:
  case 208:
  case 209:
    m_cplcorr = 
      pow(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_QED"))/
	  MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))/
	  ewcouple_.gwsq,3.);
    m_cplcorr *= (pID==115)?1./3.:1.;
    break;
  }
  m_cplcorr *=
    ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("h0_gg_fac"))/(2./3.));

  msg_Tracking()<<"Potential finite top mass correction (enabled) yields: "
		<<MODEL::s_model->ScalarConstant(std::string("h0_gg_fac"))/(2./3.)
		<<"."<<std::endl;

  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_gg_h::~MCFM_gg_h()
{
  delete [] p_p;
  delete [] p_msqv;
}


double MCFM_gg_h::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 112: gg_h_v_(p_p,p_msqv); break;
  case 113: qqb_hww_v_(p_p,p_msqv); break;
  case 114:
  case 115: qqb_hzz_v_(p_p,p_msqv); break;
  case 204: gg_hg_v_(p_p,p_msqv); break;
  case 208: gg_hwwg_v_(p_p,p_msqv); break;
  case 209: gg_hwwg_v_(p_p,p_msqv); break; 
  }
  return p_msqv[mr(i,j)];
}

void MCFM_gg_h::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);
  if (m_pID>200) corrfactor *= (*p_as)(m_mur2)/qcdcouple_.as;
  double sh((m_pID==112||m_pID==204)?
	    (p[2]+p[3]).Abs2() :
	    (p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= 
    (ATOOLS::sqr(sh-ATOOLS::sqr(masses_.hmass))+
     ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
    (ATOOLS::sqr(sh-m_mh2)+m_mh2*m_Gh2);
  for (int n(0);n<2;++n)        GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
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
}

double MCFM_gg_h::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gg_h_Getter,"MCFM_gg_h")
Virtual_ME2_Base *MCFM_gg_h_Getter::operator()(const Process_Info &pi) const
{
  msg_Out()<<"Check for process in "<<METHOD<<"."<<std::endl;
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    // check for right model and absence of b Yukawa couplings
    Flavour_Vector fl(pi.ExtractFlavours());
    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    int pID(0);
    if (pi.m_fi.m_ps.size()<1 || pi.m_fi.m_ps.size()>2) return NULL;
    ATOOLS::Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
    // higgs propagator
    if (flh!=ATOOLS::Flavour(kf_h0))                    return NULL;
    if (pi.m_fi.m_ps.size()==2 && 
	!pi.m_fi.m_ps[1].m_fl[0].Strong())              return NULL;

    if (ATOOLS::Flavour(kf_b).Yuk()>0. ||
	MODEL::s_model->Name()!=std::string("SM+EHC") ||
	!Flavour(kf_h0).IsOn()) {
      msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		 <<"   Try to initialise process gg->H(+jet) in MCFM."<<std::endl
		 <<"   Inconsistent setting with Sherpa: "<<std::endl
		 <<"Yuk(b) = "<<ATOOLS::Flavour(kf_b).Yuk()<<" (should be 0), "
		 <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM+EHC', and "
		 <<"higgs on = "<<Flavour(kf_h0).IsOn()<<"(should be 1)."
		 <<std::endl<<"   Will exit the run."<<std::endl;
      exit(1);
      return NULL;
    }
    // tau tau final state
    if ((fl.size()==4 || fl.size()==5) && 
	(fl[2]==fl[3].Bar() && fl[2].Kfcode()==15)) {
      if (ATOOLS::Flavour(kf_tau).Yuk()<=0.) {
	msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		   <<"   Setup for gg->[h->tau tau] (+jet), but tau Yukawa = 0."
		   <<std::endl;
	exit(1);
      }
      if (fl.size()==4 && pi.m_fi.m_ps.size()==1 && 
	  fl[0].IsGluon() && fl[1].IsGluon()) pID = 112;
      if (fl.size()==5 && pi.m_fi.m_ps.size()==2 && 
	  fl[0].Strong() && fl[1].Strong() &&
	  pi.m_fi.m_ps[1].m_fl[0].Strong())   pID = 204;
      // consider extra jet flavour - maybe will need more tests in future ...
    }
    // VV final states
    if ((fl.size()==6 || fl.size()==7)) {
      // check for two propagators off the Higgs decay
      if (pi.m_fi.m_ps[0].m_ps.size()!=2)               return NULL; 
      // check for fully leptonic FS
      if (!(fl[2].IsLepton() && fl[3].IsLepton() && 
	    fl[4].IsLepton() && fl[5].IsLepton()))      return NULL;
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[0].m_ps[1].m_fl[0]);
      // WW final state
      if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	  (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	if (ATOOLS::Flavour(kf_Wplus).Yuk()<=0.) {
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"   Setup for gg->[h->WW] (+jet), but W Yukawa = 0."
		     <<std::endl;
	  exit(1);
	}
	if (fl.size()==6 && pi.m_fi.m_ps.size()==1 && 
	    fl[0].IsGluon() && fl[1].IsGluon()) pID = 113;
	if (fl.size()==7 && pi.m_fi.m_ps.size()==2 && 
	    fl[0].Strong() && fl[1].Strong() &&
	    pi.m_fi.m_ps[1].m_fl[0].Strong())   pID = 208;
	// consider extra jet flavour - maybe will need more tests in future ...
      }
      // ZZ final state
      else if (fl1==Flavour(kf_Z) && fl2==Flavour(kf_Z)) {
	if (ATOOLS::Flavour(kf_Z).Yuk()<=0.) {
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"   Setup for gg->[h->ZZ] (+jet), but Z Yukawa = 0."
		     <<std::endl;
	  exit(1);
	}
	int neutrino(0);
	if (fl[2].IsUptype() && fl[4].IsUptype()) {
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"   Setup for gg->[h->ZZ] with 4 nu FS."<<std::endl
		     <<"   not implemented in MCFM."<<std::endl;
	  return NULL;
	}
	if ((fl[2].IsUptype() && fl[4].IsDowntype()) ||
	    (fl[2].IsDowntype() && fl[4].IsUptype())) neutrino=1;
	if (fl.size()==6 && pi.m_fi.m_ps.size()==1 && 
	    fl[0].IsGluon() && fl[1].IsGluon()) pID = 114+neutrino;
	if (fl.size()==7 && pi.m_fi.m_ps.size()==2 && 
	    fl[0].Strong() && fl[1].Strong() &&
	    pi.m_fi.m_ps[1].m_fl[0].Strong()) {
	  if (neutrino==1) {
	    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		       <<"   Setup for gg->[h->ZZ]+jet with 2 nu FS."<<std::endl
		       <<"   not implemented in MCFM."<<std::endl;
	    return NULL;
	  }
	  pID = 209;
	}
	// consider extra jet flavour - maybe will need more tests in future ...
      }
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
      return new MCFM_gg_h(pID,pi,fl);
    }
  }
  return NULL;
}
