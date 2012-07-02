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
    int                     m_pID,m_me;
    int                     m_j1,m_j2,m_j3,m_j4,m_k1,m_k2,m_k3,m_k4;
    double                * p_p, *p_msqv;
    MODEL::Running_AlphaS * p_as;
    double                  m_vev,m_mh2,m_Gh2,m_mZ2,m_GZ2,m_mW2,m_GW2;
    double                  m_ehcscale2,m_cplcorr,m_normcorr;

    void SelectIndices();
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
  void   spinoru_(const int & N,double *p,Complex * za,Complex * zb);
  double hggggvsqanal_(int & j1,int & j2,int & j3,int & j4);
  double haqggvsqanal_(int & j1,int & j2,int & j3,int & j4);
  double hqarbvsqanal_(int & j1,int & j2,int & j3,int & j4);
  double hqaqavsqanal_(int & j1,int & j2,int & j3,int & j4);
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
  p_as((MODEL::Running_AlphaS *)
       MODEL::s_model->GetScalarFunction(std::string("alpha_S"))),
  m_vev(MODEL::s_model->ScalarConstant(std::string("vev"))),
  m_mh2(sqr(Flavour(kf_h0).Mass())),
  m_Gh2(sqr(Flavour(kf_h0).Width())),
  m_mZ2(sqr(Flavour(kf_Z).Mass())),
  m_GZ2(sqr(Flavour(kf_Z).Width())),
  m_mW2(sqr(Flavour(kf_Wplus).Mass())),
  m_GW2(sqr(Flavour(kf_Wplus).Width())),
  m_ehcscale2(MODEL::s_model->ScalarConstant(std::string("EHC_SCALE2"))),
  m_cplcorr(sqr((*p_as)(m_ehcscale2)/(12.*m_vev) *
		MODEL::s_model->ScalarConstant(std::string("h0_gg_fac"))/
		(2./3.))),
  m_normcorr(4.*9.)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  switch (m_pID) {
  case 272:
    msg_Out()<<METHOD<<"(pID = 272), eff vertex = "
	     <<MODEL::s_model->ScalarConstant(std::string("h0_gg_fac"))<<".\n";
    m_cplcorr *= 
      sqr(Flavour(kf_tau).Yuk()/
	  MODEL::s_model->ScalarConstant(std::string("vev")));
    break;
  case 273:
  case 274:
    exit(1);
    break;
  }

  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
  SelectIndices();
}

MCFM_gg_hgg::~MCFM_gg_hgg()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_gg_hgg::SelectIndices() {
  size_t ng(0),flsize(m_flavs.size());

  for (size_t i=0;i<flsize;i++) {
    if (m_flavs[i].IsGluon()) ng++;
  }

  m_j1 = 1; m_j2 = 2;  m_j3 = flsize-1; m_j4 = flsize;
  if (ng==4) {
    m_me = 0;
    m_normcorr *= sqr(8./3.);
  }
  else if (ng==2) {
    m_me = 1;
    if (m_flavs[0].IsQuark() && m_flavs[1].IsQuark() &&
	m_flavs[0]==m_flavs[1].Bar()) {
      if (m_flavs[1].IsAnti()) { m_j1 = 1; m_j2 = 2; }
                          else { m_j1 = 2; m_j2 = 1; }
      m_normcorr *= sqr(8./3.);
    }
    else if (m_flavs[flsize-2].IsQuark() && m_flavs[flsize-1].IsQuark() &&
	     m_flavs[flsize-2]==m_flavs[flsize-1].Bar()) {
      if (m_flavs[flsize-1].IsAnti()) { m_j1 = flsize;   m_j2 = flsize-1; }
                                 else { m_j1 = flsize-1; m_j2 = flsize;   }
      m_j3 = 1; m_j4 = 2;
      m_normcorr *= sqr(8./3.);
    }
    else {
      for (size_t i=0;i<2;i++) {
	for (size_t j=0;j<2;j++) {
	  if (m_flavs[i].IsGluon() && m_flavs[flsize-2+j].IsGluon() &&
	      m_flavs[1-i].IsQuark() && m_flavs[flsize-1-j].IsQuark() &&
	      m_flavs[1-i]==m_flavs[flsize-1-j]) {
	    if (m_flavs[1-i].IsAnti()) { m_j1 = flsize-j; m_j2 = 2-i; }
	                          else { m_j1 = 2-i;      m_j2 = flsize-j; }
	    m_j3 = i+1; m_j4 = flsize-1+j;    
	    m_normcorr *= sqr(8./3.);
	  }
	}
      }
    }
  }
  else if (ng==0) {
    for (size_t i=0;i<2;i++) {
      for (size_t j=0;j<2;j++) {
	if (m_flavs[i] == m_flavs[flsize-2+j] &&
	    m_flavs[1-i]==m_flavs[flsize-1-j]) {
	  if (m_flavs[i].IsAnti()) { m_j3 = i+1;        m_j4 = flsize-1+j; }
	                      else { m_j3 = flsize-1+j; m_j4 = i+1; }
	  if (m_flavs[1-i].IsAnti()) { m_j1 = 2-i;      m_j2 = flsize-j; }
	                       else  { m_j1 = flsize-j; m_j2 = 2-i; }
	  m_me = m_flavs[i].Kfcode()!=m_flavs[1-i].Kfcode()?2:3;
	  if (m_me==3) { 
	    if (m_flavs[i].IsAnti()) {
	      m_k1 = m_j1; m_k2 = m_j4; m_k3 = m_j3; m_k4 = m_j2;
	    }
	    else {
	      m_k1 = m_j3; m_k2 = m_j2; m_k3 = m_j1; m_k4 = m_j4;
	    }
	  }
	  m_normcorr *= sqr(8./3.);
	  i = 2; j = 2; break;
	}
	else if (m_flavs[i]==m_flavs[1-i].Bar() &&
		 m_flavs[flsize-2+j]==m_flavs[flsize-1-j].Bar()) {
	  if (m_flavs[i].IsAnti()) { m_j3 = i+1;          m_j4 = 2-i; }
	                      else { m_j3 = 2-i;          m_j4 = i+1; }
	  if (m_flavs[1-i].IsAnti()) { m_j1 = flsize-1+j; m_j2 = flsize-j; }
	                       else  { m_j1 = flsize-j;   m_j2 = flsize-1+j; }
	  m_me = m_flavs[i].Kfcode()!=m_flavs[flsize-2+j].Kfcode()?2:3;
	  if (m_me==3) { 
	    if (m_flavs[i].IsAnti()) {
	      m_k1 = m_j1; m_k2 = m_j4; m_k3 = m_j3; m_k4 = m_j2;
	    }
	    else {
	      m_k1 = m_j3; m_k2 = m_j2; m_k3 = m_j1; m_k4 = m_j4;
	    }
	  }
	  m_normcorr *= sqr(8./3.);
	  i = 2; j = 2; break;
	}
      }
    }
  }
}

void MCFM_gg_hgg::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);
  corrfactor *= sqr((*p_as)(m_mur2));
  // propagator corrections
  double sh(m_pID==272?
	    (p[2]+p[3]).Abs2() :
	    (p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= 1./(sqr(sh-m_mh2)+m_mh2*m_Gh2);
  switch (m_pID) {
  case 272:
    corrfactor *= 2.*sh;
  }

  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);

  spinoru_(p.size(),p_p,zprods_.za,zprods_.zb);

  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  double res,res1,res2;
  switch(m_me) {
  case 0:
    epinv_.epinv=epinv2_.epinv2=0.0;
    res  = hggggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv_.epinv=1.0;
    res1 = hggggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv2_.epinv2=1.0;
    res2 = hggggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    break;
  case 1:
    epinv_.epinv=epinv2_.epinv2=0.0;
    res  = haqggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv_.epinv=1.0;
    res1 = haqggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv2_.epinv2=1.0;
    res2 = haqggvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    break;
  case 2:
    epinv_.epinv=epinv2_.epinv2=0.0;
    res  = hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv_.epinv=1.0;
    res1 = hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    epinv2_.epinv2=1.0;
    res2 = hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) * corrfactor;
    break;
  case 3:
    epinv_.epinv=epinv2_.epinv2=0.0;
    res  = (hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) +
	    hqarbvsqanal_(m_k1,m_k2,m_k3,m_k4) +
	    hqaqavsqanal_(m_j1,m_j2,m_j3,m_j4)) * corrfactor;
    epinv_.epinv=1.0;
    res1 = (hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) +
	    hqarbvsqanal_(m_k1,m_k2,m_k3,m_k4) +
	    hqaqavsqanal_(m_j1,m_j2,m_j3,m_j4)) * corrfactor;
    epinv2_.epinv2=1.0;
    res2 = (hqarbvsqanal_(m_j1,m_j2,m_j3,m_j4) +
	    hqarbvsqanal_(m_k1,m_k2,m_k3,m_k4) +
	    hqaqavsqanal_(m_j1,m_j2,m_j3,m_j4)) * corrfactor;
    break;
  }

  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
}

double MCFM_gg_hgg::Eps_Scheme_Factor(const Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gg_hgg_Getter,"MCFM_gg_hgg")
Virtual_ME2_Base *MCFM_gg_hgg_Getter::operator()(const Process_Info &pi) const
{
  msg_Out()<<"In "<<METHOD<<":\n";
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (MODEL::s_model->Name()!=std::string("SM+EHC") ||
      MODEL::s_model->ScalarConstant("Yukawa_b")>0. ||
      !Flavour(kf_h0).IsOn())                           return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (!(pi.m_fi.m_nloqcdtype&nlo_type::loop))           return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (!(fl[0].Strong() && fl[1].Strong()))              return NULL;
  if (pi.m_fi.m_ps.size()!=3)                           return NULL;
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
  if (fl.size()==6 && pi.m_fi.m_ps.size()==3 && 
      fl[4].Strong() && fl[5].Strong() &&
      fl[2]==fl[3].Bar() && fl[2].Kfcode()==15) {
    if (Flavour(kf_tau).Yuk()<=0.) {
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Setup for gg->[h->tau tau] (+jet), but tau Yukawa = 0."
		 <<std::endl;
      THROW(fatal_error,"Inconsistent setup.");
    }
    pID = 272;
  }
  //////////////////////////////////////////////////////////////////
  // photon photon final states 
  //////////////////////////////////////////////////////////////////
  else if (fl.size()==6 && pi.m_fi.m_ps.size()==3 && 
	   fl[4].Strong() && fl[5].Strong() &&
	   fl[2]==fl[3] && fl[2].Kfcode()==22) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   gg->[h->gamma gamma] not yet debugged.\n";
    THROW(fatal_error,"Not yet working.");
    if (fl.size()==4 && pi.m_fi.m_ps.size()==1)           pID = 210;
  }
  //////////////////////////////////////////////////////////////////
  // Z Z final states 
  //////////////////////////////////////////////////////////////////
  else if (fl.size()==8 && pi.m_fi.m_ps.size()==3 && 
	   fl[6].Strong() && fl[7].Strong() &&
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
    if (Z11.Charge()!=0 && Z21.Charge()!=0)      pID = 274;
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
  else if (fl.size()==8 && pi.m_fi.m_ps.size()==3 && 
	   fl[6].Strong() && fl[7].Strong() &&
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
	 (W21.Kfcode()==14 && W22.Kfcode()==13)))        pID = 273;
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
  return NULL;
}
