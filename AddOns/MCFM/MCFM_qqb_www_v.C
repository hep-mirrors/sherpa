#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  // README:
  // For Higgs production, choose model: MODEL = SM+EHC
  // It is important for the Higgs production to have all five flavours 
  // in the initial state, but the Yukawa coupling of the b must be
  // switched off:  YUKAWA_B = 0.  
  // Also, MCFM acts in the limit of mt->infinity,
  // thus either a correction term must be introduced, or we have
  // to set: FINITE_TOP_MASS = 0 
  // For direct W-pair production, there's a contribution from b
  // initial states, through a top in the t-channel.  Switch it off
  // by setting   ACTIVE[6]=0 and/or MASSIVE[5]=1;



  class MCFM_qqb_www_v: public PHASIC::Virtual_ME2_Base {
  private:
    bool    m_isW, m_isneutrino;
    double *p_p, *p_msqv;
    double  m_mh2,m_Gh2,m_vev,m_aqed,m_sin2tw,m_asmh,m_ewcorr;
  public:
    MCFM_qqb_www_v(const PHASIC::Process_Info& pi,
		   const ATOOLS::Flavour_Vector& flavs,
		   const bool & isW,const bool & isneutrino);
    ~MCFM_qqb_www_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_wh_ww_v_(double *p,double *msqv); 
  void qqb_zh_ww_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_www_v::MCFM_qqb_www_v(const Process_Info& pi,
			       const Flavour_Vector& flavs,
			       const bool & isW,const bool & isneutrino):
  Virtual_ME2_Base(pi,flavs), m_isW(isW), m_isneutrino(isneutrino),
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_vev(MODEL::s_model->ScalarConstant(std::string("vev"))),
  m_aqed(MODEL::s_model->ScalarFunction(std::string("alpha_QED"))),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),
  m_asmh(MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mh2)),
  m_ewcorr(pow((4.*M_PI*m_aqed/m_sin2tw)/ewcouple_.gwsq,6.) *
	   (m_isW?1.:ATOOLS::sqr((m_sin2tw/(1.-m_sin2tw))/
				 (ewcouple_.xw/(1.-ewcouple_.xw)))) *
	   (m_isW?pow(ATOOLS::Flavour(kf_Wplus).Yuk()/masses_.wmass,2.):
	    pow(ATOOLS::Flavour(kf_Z).Yuk()/
		(masses_.wmass/sqrt(1.-ewcouple_.xw)),2.)))
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_www_v::~MCFM_qqb_www_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_www_v::Calc(const Vec4D_Vector &p)
{
  double sf(4.0*9.0/qcdcouple_.ason2pi*(m_isneutrino?1./3.:1.));
  double cplfactor(m_ewcorr);

  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) {
    if (n<6)             GetMom(p_p,n+2,p[n]);
    else if (n>5)        GetMom(p_p,n-4,p[n]);
  }
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  double s2345((p[2]+p[3]+p[4]+p[5]).Abs2());
  double propfactor = 
      (ATOOLS::sqr(s2345-ATOOLS::sqr(masses_.hmass))+
       ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
      (ATOOLS::sqr(s2345-m_mh2)+m_mh2*m_Gh2);
  //msg_Out()<<METHOD<<"("<<m_isW<<"/"<<m_isneutrino<<") --> "
  //	   <<"M_2345 = "<<sqrt(s2345)<<", "
  //	   <<"M_67 = "<<sqrt((p[6]+p[7]).Abs2())<<", "
  //	   <<"ewcorr = "<<m_ewcorr<<" & "
  //	   <<"xw = "<<ewcouple_.xw<<" vs "
  //	   <<MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))
  //	   <<"."<<std::endl;

  epinv_.epinv=epinv2_.epinv2=0.0;
  if (m_isW) qqb_wh_ww_v_(p_p,p_msqv);
        else qqb_zh_ww_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  if (m_isW) qqb_wh_ww_v_(p_p,p_msqv);
        else qqb_zh_ww_v_(p_p,p_msqv);

  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  if (m_isW) qqb_wh_ww_v_(p_p,p_msqv);
        else qqb_zh_ww_v_(p_p,p_msqv);

  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res         * cplfactor * propfactor;
  m_res.IR()     = (res1-res)  * cplfactor * propfactor;
  m_res.IR2()    = (res2-res1) * cplfactor * propfactor;
}

double MCFM_qqb_www_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_www_v_Getter,"MCFM_qqb_www_v")
Virtual_ME2_Base *MCFM_qqb_www_v_Getter::operator()(const Process_Info &pi) const
{
  //msg_Out()<<METHOD<<"===================="<<std::endl;
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=8) return NULL;
    int pID(0), swapped(0);
    //msg_Out()<<"   "<<fl.size()<<" "<<pi.m_fi.m_ps.size()<<"."<<std::endl;
    if (fl[0].IsQuark() && fl[1].IsQuark() &&
	fl[2].IsLepton() && fl[3].IsLepton() && 
	fl[4].IsLepton() && fl[5].IsLepton() &&
	fl[6].IsLepton() && fl[7].IsLepton() &&
	pi.m_fi.m_ps.size()==2) {
      ATOOLS::Flavour flV = pi.m_fi.m_ps[0].m_fl[0];
      ATOOLS::Flavour flh = pi.m_fi.m_ps[1].m_fl[0];
      if (flV==ATOOLS::Flavour(kf_h0) && 
	  (flh==ATOOLS::Flavour(kf_Wplus) || 
	   flh==ATOOLS::Flavour(kf_Wplus).Bar() ||
	   flh==ATOOLS::Flavour(kf_Z))) {
	flV     = pi.m_fi.m_ps[1].m_fl[0];
	flh     = pi.m_fi.m_ps[0].m_fl[0];
	swapped = 1;
      }
      //msg_Out()<<"  check: "<<flh<<" "<<flV<<" "<<swapped<<"."<<std::endl;
      if (flh==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	  pi.m_fi.m_ps[1-swapped].m_ps.size()==2) {
	ATOOLS::Flavour fl1(pi.m_fi.m_ps[1-swapped].m_ps[0].m_fl[0]);
	ATOOLS::Flavour fl2(pi.m_fi.m_ps[1-swapped].m_ps[1].m_fl[0]);
	//msg_Out()<<"  "<<pi.m_fi.m_ps[1-swapped].m_ps.size()<<"  "
	//	 <<fl1<<" & "<<fl2<<"  --> "
	//	 <<pi.m_fi.m_ps[swapped].m_ps[0].m_fl[0]<<" & "
	//	 <<pi.m_fi.m_ps[swapped].m_ps[1].m_fl[0]<<std::endl;
	if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	    (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	  if (flV==ATOOLS::Flavour(kf_Wplus))                       pID = 92;
	  else if (flV==ATOOLS::Flavour(kf_Wplus).Bar())            pID = 97;
	  else if (flV==ATOOLS::Flavour(kf_Z) &&
		   pi.m_fi.m_ps[swapped].m_ps[0].m_fl[0]==
		   pi.m_fi.m_ps[swapped].m_ps[1].m_fl[0].Bar()) {
	    if (pi.m_fi.m_ps[swapped].m_ps[0].m_fl[0].IsUptype())   pID = 107;
	    if (pi.m_fi.m_ps[swapped].m_ps[0].m_fl[0].IsDowntype()) pID = 106;
	  }
	}
      }
    }
    if (pID>0) {
      removebr_.removebr=false;
      zerowidth_.zerowidth=true;
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      nproc_.nproc=pID;
      chooser_();
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_qqb_www_v(pi,fl,pID==92||pID==97,pID==107);
    }
  }
  return NULL;
}
