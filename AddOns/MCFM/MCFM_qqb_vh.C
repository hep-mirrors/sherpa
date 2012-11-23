#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  class MCFM_qqb_vh: public PHASIC::Virtual_ME2_Base {
  private:
    int     m_pID;
    double *p_p, *p_msqv;
    double  m_mh2,m_Gh2,m_vev,m_aqed,m_sin2tw;
    double  m_normcorr,m_ewcorr;

    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_qqb_vh(const int & pID,const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_vh();
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

MCFM_qqb_vh::MCFM_qqb_vh(const int & pID,const Process_Info& pi,
			 const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), m_pID(pID), 
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_vev(MODEL::s_model->ScalarConstant(std::string("vev"))),
  m_aqed(MODEL::s_model->ScalarFunction(std::string("alpha_QED"))),
  m_sin2tw(MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))),
  m_normcorr(4.0*9.0/qcdcouple_.ason2pi/((m_pID==107)?3.:1.)),
  m_ewcorr(pow((4.*M_PI*m_aqed/m_sin2tw)/ewcouple_.gwsq,6.))
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM.");
  if (m_pID==92 || m_pID==97) {
    m_ewcorr *= sqr(Flavour(kf_Wplus).Yuk()/masses_.wmass);
  }
  else if (m_pID==106 || m_pID==107) {
    m_ewcorr *= sqr(Flavour(kf_Z).Yuk()/(masses_.wmass/sqrt(1.-ewcouple_.xw)));
    m_ewcorr *= sqr((m_sin2tw/(1.-m_sin2tw))/(ewcouple_.xw/(1.-ewcouple_.xw))); 
  }
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_vh::~MCFM_qqb_vh()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_vh::CallMCFM(const int & i,const int & j) {
  switch (m_pID) {
  case 92:
  case 97:  qqb_wh_ww_v_(p_p,p_msqv); break;
  case 106:
  case 107: qqb_zh_ww_v_(p_p,p_msqv); break;
  }
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_vh::Calc(const Vec4D_Vector &p)
{
  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) {
    if (n<6)             GetMom(p_p,n+2,p[n]);        
    else if (n>5)        GetMom(p_p,n-4,p[n]);
  }
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  double corrfactor(m_ewcorr*m_normcorr);
  double s2345((p[2]+p[3]+p[4]+p[5]).Abs2());
  corrfactor *= 
      (ATOOLS::sqr(s2345-ATOOLS::sqr(masses_.hmass))+
       ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
      (ATOOLS::sqr(s2345-m_mh2)+m_mh2*m_Gh2);

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

double MCFM_qqb_vh::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_vh_Getter,"MCFM_qqb_vh")
Virtual_ME2_Base *MCFM_qqb_vh_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (MODEL::s_model->Name()!=std::string("SM")) return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=8) return NULL;
    int pID(0);
    if (fl[0].IsQuark() && fl[1].IsQuark() &&
	fl[2].IsLepton() && fl[3].IsLepton() && 
	fl[4].IsLepton() && fl[5].IsLepton() &&
	fl[6].IsLepton() && fl[7].IsLepton() &&
	pi.m_fi.m_ps.size()==2) {
      ATOOLS::Flavour flh = pi.m_fi.m_ps[0].m_fl[0];
      ATOOLS::Flavour flV = pi.m_fi.m_ps[1].m_fl[0];
      msg_Out()<<"  check: "<<flh<<" "<<flV<<"."<<std::endl;
      if (flh==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	  pi.m_fi.m_ps[1].m_ps.size()==2) {
	ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_ps[0].m_fl[0]);
	ATOOLS::Flavour fl2(pi.m_fi.m_ps[0].m_ps[1].m_fl[0]);
	if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	    (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	  if (Flavour(kf_b).Yuk()>0. ||
	      MODEL::s_model->Name()!=std::string("SM")) {
	    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		       <<"   Try to initialise process qqb->VH in MCFM.\n"
		       <<"   Inconsistent setting with Sherpa: \n"
		       <<"YUKAWA[5] = "<<Flavour(kf_b).Yuk()
		       <<" (should be 0), and "
		       <<"MODEL = "<<MODEL::s_model->Name()
		       <<"(should be 'SM'.\n"
		       <<"   Will exit the run."<<std::endl;
	    THROW(not_implemented,"Incompatible setting with MCFM.");
	    return NULL;
	  }
	  if (flV==ATOOLS::Flavour(kf_Wplus))                       pID = 92;
	  else if (flV==ATOOLS::Flavour(kf_Wplus).Bar())            pID = 97;
	  else if (flV==ATOOLS::Flavour(kf_Z) &&
		   pi.m_fi.m_ps[1].m_ps[0].m_fl[0]==
		   pi.m_fi.m_ps[1].m_ps[1].m_fl[0].Bar()) {
	    if (pi.m_fi.m_ps[1].m_ps[0].m_fl[0].IsUptype())   pID = 107;
	    if (pi.m_fi.m_ps[1].m_ps[0].m_fl[0].IsDowntype()) pID = 106;
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
      return new MCFM_qqb_vh(pID,pi,fl);
    }
  }
  return NULL;
}
