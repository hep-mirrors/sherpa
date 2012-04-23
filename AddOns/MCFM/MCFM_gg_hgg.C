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


  class MCFM_gg_hgg: public PHASIC::Virtual_ME2_Base {
  private:
    int                     m_pID;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    double                  m_mh2,m_Gh2,m_ehcscale2,m_cplcorr,m_normcorr;
  public:
    MCFM_gg_hgg(const int & pID,
		const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_gg_hgg();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void spinoru_(const int & N,double *p,
		Complex * za,Complex * zb);
  double hggggvsqanal_(int & j1,int & j2,int & j3,int & j4);
  void gg_hgg_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;

MCFM_gg_hgg::MCFM_gg_hgg(const int & pID,const Process_Info& pi,
			 const Flavour_Vector& flavs) :
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  p_as((Running_AlphaS *)s_model->GetScalarFunction(string("alpha_S"))),
  m_mh2(sqr(Flavour(kf_h0).Mass())),
  m_Gh2(sqr(Flavour(kf_h0).Width())),
  m_ehcscale2(s_model->ScalarConstant(string("EHC_SCALE2"))),
  m_cplcorr(1.),m_normcorr(1.)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  m_cplcorr = 
    sqr(((*p_as)(m_ehcscale2) *
	 s_model->ScalarConstant(string("h0_gg_fac")))/
	(3.*M_PI*s_model->ScalarConstant(string("vev")))) *
    sqr(Flavour(kf_tau).Yuk()/s_model->ScalarConstant(string("vev")));

  msg_Out()<<METHOD<<": cplcorr = "<<m_cplcorr<<".\n";
  
  
  switch (m_pID) {
  case 272:
    break;
  }
  msg_Tracking()<<"Potential finite top mass correction (enabled) yields: "
		<<(s_model->ScalarConstant(string("h0_gg_fac")))<<".\n";

  p_p      = new double[4*MCFM_NMX];
  p_msqv   = new double[sqr(2*MCFM_NF+1)];
  m_drmode = m_mode=1;
}

MCFM_gg_hgg::~MCFM_gg_hgg()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_gg_hgg::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);
  corrfactor *= pow(3.*(*p_as)(m_mur2),2)*sqr(4.*M_PI)/2.;  
  // factor in the colour factors plus the symmetrisation
  double sh((m_pID==112||m_pID==204||m_pID==272)?
	    (p[2]+p[3]).Abs2() :
	    (p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= sh/(sqr(sh-m_mh2)+m_mh2*m_Gh2);
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);

  spinoru_(p.size(),p_p,zprods_.za,zprods_.zb);
  //gg_hgg_v_(p_p,p_msqv);
  //std::cout<<std::endl;
  //exit(1);

  long int i(m_flavs[0]), j(m_flavs[1]);
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  double res,res1,res2;
  int j1,j2,j3,j4;
  if (m_flavs[0].IsGluon() && 
      m_flavs[1].IsGluon()) {
    if (m_flavs[4].IsGluon() && 
	m_flavs[5].IsGluon()) {
      j1 = 1;
      j2 = 2;
      j3 = 5;
      j4 = 6;
      epinv_.epinv=epinv2_.epinv2=0.0;
      res  = hggggvsqanal_(j1,j2,j3,j4) * corrfactor;
      epinv_.epinv=1.0;
      res1 = hggggvsqanal_(j1,j2,j3,j4) * corrfactor;
      epinv2_.epinv2=1.0;
      res2 = hggggvsqanal_(j1,j2,j3,j4) * corrfactor;
    }
  }

  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);

  msg_Out()<<METHOD<<" yields "<<m_res.Finite()
  	   <<" + 1/eps * "<<m_res.IR()
  	   <<" + 1/eps^2 * "<<m_res.IR2()
	   <<" for mb = "<<masses_.mb
	   <<" and "<<nflav_.nflav<<" active flavours"
	   <<" in "<<scheme_.scheme<<", corr = "<<corrfactor<<".\n";
  msg_Out()<<"sqrt(sh) = "<<sqrt(sh)<<", mh = "<<sqrt(m_mh2)<<", "
	   <<"1/prop = "<<(sqr(sh-m_mh2)+m_mh2*m_Gh2)<<"\n"
	   <<"coupling = "<<m_cplcorr<<" from "
	   <<"yuk = "<<Flavour(kf_tau).Yuk()<<", "
	   <<"v = "<<s_model->ScalarConstant(string("vev"))<<", "
	   <<"eff = "<<s_model->ScalarConstant(string("h0_gg_fac"))<<" & "
	   <<"as = "<<(*p_as)(m_ehcscale2)<<" & "<<(*p_as)(m_mur2)<<".\n"
	   <<"cplcorr = "<<m_cplcorr<<".\n";
}

double MCFM_gg_hgg::Eps_Scheme_Factor(const Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gg_hgg_Getter,"MCFM_gg_hgg")
Virtual_ME2_Base *MCFM_gg_hgg_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (s_model->Name()!=string("SM+EHC") ||
      s_model->ScalarConstant("Yukawa_b")>0. ||
      !Flavour(kf_h0).IsOn())                           return NULL;
  if (pi.m_oew>2)                                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    msg_Out()<<METHOD<<" checks numbers: "
    	     <<fl.size()<<" external particles, "
    	     <<pi.m_fi.m_ps.size()<<" props";
    if (pi.m_fi.m_ps.size()==0) msg_Out()<<".\n";
    else msg_Out()<<" and "
    		  <<pi.m_fi.m_ps[0].m_ps.size()<<" final state particles.\n";
    //msg_Out()<<"further check: "
    //	     <<pi.m_fi.m_ps[1].m_fl[0]<<"\n";
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    if (!(fl.size()>=6 &&
	  pi.m_fi.m_ps.size()>=3 &&
	  pi.m_fi.m_ps[0].m_ps.size()>1))               return NULL;
    Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
    if (!flh==Flavour(kf_h0))                   return NULL;
    int pID(0);
    // tau tau final states
    if (fl[2]==fl[3].Bar() && fl[2].Kfcode()==15) {
      if (Flavour(kf_tau).Yuk()<=0.) {
	msg_Error()<<"Error in "<<METHOD<<":"<<endl
		   <<"   Setup for gg->[h->tau tau] (+jet), but tau Yukawa = 0."
		   <<endl;
        THROW(fatal_error,"Inconsistent setup.");
      }
      if (fl.size()==6 && pi.m_fi.m_ps.size()==3 &&
	  pi.m_fi.m_ps[1].m_fl[0].Strong()&&
	  pi.m_fi.m_ps[2].m_fl[0].Strong()) {
	pID = 272;
      }
    }
    PRINT_VAR(pID);
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
      return new MCFM_gg_hgg(pID,pi,fl);
    }
  }
  return NULL;
}
