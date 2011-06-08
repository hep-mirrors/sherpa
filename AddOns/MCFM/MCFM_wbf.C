#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MCFM_Wrapper.H"

namespace MCFM {


  class MCFM_wbf: public PHASIC::Virtual_ME2_Base {
  private:
    int                     m_pID;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    double                  m_mh2,m_Gh2,m_cplcorr,m_normcorr;
    
    
    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_wbf(const int & pID,
	     const PHASIC::Process_Info& pi,
	     const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_wbf();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };
  
}// end of namespace MCFM

extern "C" { 
  void vv_hqq_v_(double *p,double *msqv); 
  void vv_hww_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_wbf::MCFM_wbf(const int & pID,const Process_Info& pi,
		   const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_cplcorr(ewcouple_.vevsq/
	    ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev")))),
  m_normcorr(4.*9./qcdcouple_.ason2pi)
{
  switch (m_pID) {
  case 212:
  case 217:
    m_cplcorr *= 
      ATOOLS::sqr(ATOOLS::Flavour(kf_tau).Yuk())/(3.*masses_.mbsq) *
      2.*4.*ATOOLS::sqr(masses_.wmass)/ewcouple_.gwsq/
      ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev")));
    break;
  case 213:
    m_cplcorr = 
      pow(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_QED"))/
	  MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))/
	  ewcouple_.gwsq,3.);
    break;
  }

  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_wbf::~MCFM_wbf()
{
  delete [] p_p;
  delete [] p_msqv;
}


double MCFM_wbf::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 212: vv_hqq_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_wbf::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);
  double sh((m_pID==212||m_pID==217)?
	    (p[2]+p[3]).Abs2() :
	    (p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= 
    (ATOOLS::sqr(sh-ATOOLS::sqr(masses_.hmass))+
     ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
    (ATOOLS::sqr(sh-m_mh2)+m_mh2*m_Gh2);
  corrfactor *= sh/(sh-2.*masses_.mbsq);

  for (int n(0);n<2;++n)        GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; corrfactor *= 8./3.; }
  if (j==21) { j=0; corrfactor *= 8./3.; }

  for (int k=0;k<6;k++) msg_Out()<<"  "<<m_flavs[k];msg_Out()<<".\n";
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

  msg_Out()<<METHOD
	   <<"IR2 = "<<m_res.IR2()<<", "
	   <<"IR1 = "<<m_res.IR()<<", "
	   <<"Fin = "<<m_res.Finite()<<".\n";
}

double MCFM_wbf::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_wbf_Getter,"MCFM_wbf")
Virtual_ME2_Base *MCFM_wbf_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (!ATOOLS::Flavour(kf_h0).IsOn())                   return NULL;
  if (pi.m_oew<3)                                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    // check for right model and absence of b Yukawa couplings
    Flavour_Vector fl(pi.ExtractFlavours());
    msg_Out()<<"  "<<fl.size()<<" flavours, "<<fl[2]<<" = "<<fl[3]<<", "
	     <<pi.m_fi.m_ps.size()<<" props.\n";

    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    int pID(0);
    if (pi.m_fi.m_ps.size()!=3)                         return NULL;
    ATOOLS::Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
    msg_Out()<<"  Higgs candidate = "<<flh<<"."<<std::endl;
    // higgs propagator
    if (!flh==ATOOLS::Flavour(kf_h0))                   return NULL;
    if (pi.m_fi.m_ps.size()==2 && 
	!pi.m_fi.m_ps[1].m_fl[0].Strong())              return NULL;

    if (MODEL::s_model->ScalarConstant("Yukawa_b")>0. ||
	MODEL::s_model->Name()!=std::string("SM") ||
	!Flavour(kf_h0).IsOn()) {
      msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		 <<"   Try to initialise process WBF->H in MCFM.\n"
		 <<"   Inconsistent setting with Sherpa: \n"
		 <<"YUKAWA_B = "<<MODEL::s_model->ScalarConstant("Yukawa_b")<<" (should be 0), "
		 <<"MODEL = "<<MODEL::s_model->Name()
		 <<"(should be 'SM'), and "
		 <<"ACTIVE[25] = "<<Flavour(kf_h0).IsOn()<<"(should be 1)."
		 <<std::endl<<"   Will exit the run.\n";
      exit(1);
      return NULL;
    }
    // tau tau final state
    if ((fl.size()>=6) && 
	(fl[2]==fl[3].Bar() && fl[2].Kfcode()==15)) {
      if (ATOOLS::Flavour(kf_tau).Yuk()<=0.) {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   Setup for WBF->[h->tau tau], but tau Yukawa = 0.\n";
	exit(1);
      }
      if (fl.size()==6 && pi.m_fi.m_ps.size()==3 && 
	  fl[0].IsQuark() && fl[1].IsQuark() && 
	  fl[4].IsQuark() && fl[5].IsQuark())   pID = 212;
      else if (fl.size()==7 && pi.m_fi.m_ps.size()==4 && 
	       fl[0].Strong() && fl[1].Strong() && 
	       fl[4].Strong() && fl[5].Strong() &&
	       fl[6].Strong())                  pID = 217;
    }
    if (fl.size()==8) {
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
	if (pi.m_fi.m_ps.size()==3 && 
	    fl[0].IsQuark() && fl[1].IsQuark() && 
	    fl[6].IsQuark() && fl[7].IsQuark())   pID = 213;
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
      return new MCFM_wbf(pID,pi,fl);
    }
  }
  return NULL;
}
