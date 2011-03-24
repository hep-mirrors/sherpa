#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  class MCFM_qqb_vv: public PHASIC::Virtual_ME2_Base {
  private:
    int     m_pID;
    bool    m_swapped;
    double  m_aqed,m_cplcorr, m_normcorr;
    double *p_p, *p_msqv;

    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_qqb_vv(const int & pID,const bool & swapped, 
		const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_vv();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_ww_v_(double *p,double *msqv); 
  void qqb_wz_v_(double *p,double *msqv); 
  void qqb_zz_v_(double *p,double *msqv); 
  void qqb_wgam_v_(double *p,double *msqv); 
  void qqb_zgam_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_vv::MCFM_qqb_vv(const int & pID,const bool & swapped,
			 const PHASIC::Process_Info& pi,
			 const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), m_pID(pID), m_swapped(swapped),
  m_aqed(MODEL::s_model->ScalarFunction(std::string("alpha_QED"))),
  m_cplcorr(ATOOLS::sqr(4.*M_PI*m_aqed/ewcouple_.esq)),
  m_normcorr(4.*9./qcdcouple_.ason2pi)
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
  // surprisingly, no summation over nus in WZ production.
  if (m_pID==82 || m_pID==87) m_normcorr /= 3.;
}

MCFM_qqb_vv::~MCFM_qqb_vv()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_vv::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 12: 
  case 17: qqb_wgam_v_(p_p,p_msqv); break;
  case 48: 
  case 49: qqb_zgam_v_(p_p,p_msqv); break;
  case 61: qqb_ww_v_(p_p,p_msqv); break;
  case 71:
  case 72:
  case 76:
  case 77: qqb_wz_v_(p_p,p_msqv); break;
  case 81:
  case 82: 
  case 86: 
  case 87: qqb_zz_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_vv::Calc(const Vec4D_Vector &p)
{
  double corrfactor(m_cplcorr*m_normcorr);

  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  if (m_pID==81 || m_pID==82) {
    GetMom(p_p,2,p[2]); 
    GetMom(p_p,3,p[4]); 
    GetMom(p_p,4,p[3]); 
    GetMom(p_p,5,p[5]); 
  }
  else if (m_pID==71 || m_pID==72 || m_pID==76 || m_pID==77) {
    GetMom(p_p,4,p[2]); 
    GetMom(p_p,5,p[3]); 
    GetMom(p_p,2,p[4]); 
    GetMom(p_p,3,p[5]); 
  }
  else if (m_pID==49) {
    corrfactor*=1./3.;
    if (m_swapped) {
      GetMom(p_p,2,p[3]);
      GetMom(p_p,3,p[4]);
      GetMom(p_p,4,p[2]);
    }
    else {
      for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
    }
  }
  else if (m_pID==48 || m_pID ==12 || m_pID==17) {
    GetMom(p_p,2,p[3]);
    GetMom(p_p,3,p[4]);
    GetMom(p_p,4,p[2]);
  }
  else {
    for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  }
  //msg_Out()<<"s24 = "<<sqrt((p[2]+p[4]).Abs2())<<", "
  //	   <<"s35 = "<<sqrt((p[3]+p[5]).Abs2())<<", "
  //	   <<"s23 = "<<sqrt((p[2]+p[3]).Abs2())<<", "
  //	   <<"s45 = "<<sqrt((p[4]+p[5]).Abs2())<<"; "
  //	   <<"corr = "<<m_cplcorr<<" * "<<m_normcorr
  //	   <<" = "<<corrfactor<<"."<<std::endl;

  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  epinv_.epinv=epinv2_.epinv2=0.0;
  double res(CallMCFM(i,j)  * corrfactor);
  epinv_.epinv=1.0;
  double res1(CallMCFM(i,j) * corrfactor);
  epinv2_.epinv2=1.0;
  double res2(CallMCFM(i,j) * corrfactor);

  msg_Out()<<"   --> "<<res<<" "<<res1<<" "<<res2<<"."<<std::endl;
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
}

double MCFM_qqb_vv::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_vv_Getter,"MCFM_qqb_vv")
Virtual_ME2_Base *MCFM_qqb_vv_Getter::operator()(const Process_Info &pi) const
{
  msg_Out()<<"Check for process in "<<METHOD<<"."<<std::endl;
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    if (fl.size()!=6 && fl.size()!=5)                   return NULL;
    // check for fully leptonic FS
    if (!(fl[0].IsQuark() && fl[1].IsQuark()))          return NULL;
    int pID(0);
    bool swapped(false);
    if (pi.m_fi.m_ps.size()==2) {
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
      if ((fl1.Kfcode()==23 || fl1.Kfcode()==24) &&
	  (fl2.Kfcode()==23 || fl2.Kfcode()==24)) {
	// check for right model and absence of b Yukawa couplings
	if ((ATOOLS::Flavour(kf_b).Yuk()>0. && fl1==fl2.Bar()) ||
	    MODEL::s_model->Name()!=std::string("SM") ||
	    (Flavour(kf_t).IsOn() && fl1.Kfcode()==24 && fl2.Kfcode()==24)) {
	  msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		     <<"   Try to initialise process qqb->VV in MCFM.\n"
		     <<"   Inconsistent setting with Sherpa: \n"
		     <<"Yuk(b) = "<<ATOOLS::Flavour(kf_b).Yuk()
		     <<" (should be 0), "
		     <<"model = "<<MODEL::s_model->Name()
		     <<"(should be 'SM', and "
		     <<"top on = "<<Flavour(kf_t).IsOn()
		     <<"(should be 1 for WW).\n"
		     <<"   Will exit the run."<<std::endl;
	  exit(1);
	  return NULL;
	}
      }
      // WW final state
      if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	  (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	if (fl[2].IsLepton() && fl[3].IsLepton() && 
	    fl[4].IsLepton() && fl[5].IsLepton()) {
	  pID = 61;
          zerowidth_.zerowidth=true;
	}
      }
      // ZZ final state
      else if (fl1==Flavour(kf_Z) && fl2==Flavour(kf_Z)) {
	int neutrino(0);
	if (fl[2].IsLepton() && fl[3].IsLepton() && 
	    fl[4].IsLepton() && fl[5].IsLepton()) {
	  if ((fl[2].IsUptype() && fl[4].IsDowntype()) ||
	      (fl[2].IsDowntype() && fl[4].IsUptype())) neutrino=1;
	  pID = 86+neutrino;
          zerowidth_.zerowidth=true;
	}
      }
      // W(+)Z final state
      else if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Z)) ||
	       (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Z))) {
	if (fl[2].IsDowntype()) {
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"   Z->l+l- not yet possible due to Z-gamma interference."
		     <<std::endl<<"   Will exit the run."<<std::endl;
	  exit(1);
	}
	zerowidth_.zerowidth=true;
	int neutrino(0);
	if (fl[2].IsLepton() && fl[3].IsLepton() && 
	    fl[4].IsLepton() && fl[5].IsLepton()) {
	  if ((fl[2].IsUptype() && fl[4].IsDowntype()) ||
	      (fl[2].IsDowntype() && fl[4].IsUptype())) neutrino=1;
	  pID = 71+neutrino;
	}
      }
      // W(-)Z final state
      else if ((fl1==Flavour(kf_Wplus).Bar() && fl2==Flavour(kf_Z)) ||
	       (fl2==Flavour(kf_Wplus).Bar() && fl1==Flavour(kf_Z))) {
	if (fl[2].IsDowntype()) {
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"   Z->l+l- not yet possible due to Z-gamma interference."
		     <<std::endl<<"   Will exit the run."<<std::endl;
	  exit(1);
	}
	zerowidth_.zerowidth=true;
	int neutrino(0);
	if (fl[2].IsLepton() && fl[3].IsLepton() && 
	    fl[4].IsLepton() && fl[5].IsLepton()) {
	  if ((fl[2].IsUptype() && fl[4].IsDowntype()) ||
	      (fl[2].IsDowntype() && fl[4].IsUptype())) neutrino=1;
	  pID = 76+neutrino;
	}
      }
      // Z (-> nu + nubar) + gamma 
      else if (fl1==Flavour(kf_Z) && fl2==Flavour(kf_photon)) {
	if (fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
	  if (MODEL::s_model->Name()!=std::string("SM")) {
	    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		       <<"   Try to initialise process qqb->Vgamma in MCFM."
		       <<std::endl
		       <<"   Inconsistent setting with Sherpa: "<<std::endl
		       <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM'."
		       <<std::endl<<"   Will exit the run."<<std::endl;
	    exit(1);
	    return NULL;
	  }
	  if (fl[2].IsUptype()){ 
	    pID = 49;
            zerowidth_.zerowidth=true;
	    swapped=false;
          }
        }
      }
    }
    if (pi.m_fi.m_ps.size()==3) {
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
      ATOOLS::Flavour fl3(pi.m_fi.m_ps[2].m_fl[0]);
      //msg_Out()<<"   check props: "<<fl1<<" & "<<fl2<<" & "<<fl3<<"."<<std::endl;
      if (fl1==Flavour(kf_photon)) {
      // Z + gamma final state
	if (fl[3].IsLepton() && fl[4]==fl[3].Bar()) {
	  if (MODEL::s_model->Name()!=std::string("SM")) {
	    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		       <<"   Try to initialise process qqb->Vgamma in MCFM."
		       <<std::endl
		       <<"   Inconsistent setting with Sherpa: "<<std::endl
		       <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM'."
		       <<std::endl<<"   Will exit the run."<<std::endl;
	    exit(1);
	    return NULL;
	  }
	  if (fl[3].IsUptype()){ 
	    pID = 49;
            zerowidth_.zerowidth=true;
	    swapped=true;
          }
          else {
	    pID = 48;
            zerowidth_.zerowidth=false;
	    swapped=true;
          }
        }
      // W + gamma final state
	if (fl[3].IsLepton() && fl[4].IsLepton() && fl[4]!=fl[3].Bar()) {
	  if (MODEL::s_model->Name()!=std::string("SM")) {
	    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		       <<"   Try to initialise process qqb->Vgamma in MCFM."
		       <<std::endl
		       <<"   Inconsistent setting with Sherpa: "<<std::endl
		       <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM'."
		       <<std::endl<<"   Will exit the run."<<std::endl;
	    exit(1);
	    return NULL;
	  }
	  if (fl[3].IsUptype() && fl[4].IsDowntype()){
	    pID = 12; 
            zerowidth_.zerowidth=false;
	    swapped=true;
	  }
	  else if (fl[3].IsDowntype() && fl[4].IsUptype()) {
	    pID = 17; 
            zerowidth_.zerowidth=false;
	    swapped=true;
          }
        } 
      }
    }
    else {
      // Only allow for ZZ final state - try to catch the photons
      if (fl[2].IsLepton() && fl[3].IsLepton() && 
	  fl[4].IsLepton() && fl[5].IsLepton() &&
	  fl[4]==fl[2].Bar() && fl[5]==fl[3].Bar()) {
	int neutrino(0);
	if ((fl[2].IsUptype() && fl[3].IsDowntype()) ||
	    (fl[2].IsDowntype() && fl[3].IsUptype())) neutrino=1;
	pID = 81+neutrino;
	msg_Error()<<"Error in "<<METHOD<<"(pID = "<<pID<<"):"<<std::endl
		   <<"   Continuum production with gamma interference tricky. "
		   <<"Single resonant diagrams in Sherpa, but not in MCFM."
		   <<std::endl<<"   Will abort the run."<<std::endl;
	exit(1);
      }
    } 
    if (pID!=0) {
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      nproc_.nproc=pID;
      chooser_();
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_qqb_vv(pID,swapped,pi,fl);
    }
  }
  return NULL;
}
