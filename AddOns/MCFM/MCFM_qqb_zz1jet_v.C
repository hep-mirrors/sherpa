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

  class MCFM_qqb_zz1jet_v: public PHASIC::Virtual_ME2_Base {
  private:
    bool    m_ishiggs;
    double *p_p, *p_msqv;
    double  m_mh2,m_Gh2,m_asmh,m_ewcorr;
  public:
    MCFM_qqb_zz1jet_v(const PHASIC::Process_Info& pi,
		  const ATOOLS::Flavour_Vector& flavs,
		  const bool & ishiggs);
    ~MCFM_qqb_zz1jet_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void gg_hzzg_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_zz1jet_v::MCFM_qqb_zz1jet_v(const Process_Info& pi,
			     const Flavour_Vector& flavs,
			     const bool & ishiggs) :
  Virtual_ME2_Base(pi,flavs), m_ishiggs(ishiggs), 
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_asmh(MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mh2)),
  m_ewcorr(ewcouple_.vevsq/
	   ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev"))) *
	   pow(4.*M_PI*MODEL::s_model->ScalarFunction(std::string("alpha_QED"))/
	       MODEL::s_model->ScalarConstant(std::string("sin2_thetaW"))/
	       ewcouple_.gwsq,3.))
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_zz1jet_v::~MCFM_qqb_zz1jet_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_zz1jet_v::Calc(const Vec4D_Vector &p)
{
  double sf(4.0*9.0/qcdcouple_.ason2pi); //*(m_isneutrino?1./3.:1.));
  double cplfactor(1.),propfactor(1.);
  if (m_ishiggs) {
    double assherpa(MODEL::s_model->ScalarFunction(std::string("alpha_S"),
						   m_mur2));
    cplfactor   = m_ewcorr * ATOOLS::sqr(m_asmh/qcdcouple_.as) *
      assherpa/qcdcouple_.as;
    double s12((p[0]+p[1]).Abs2());
    propfactor = 
      (ATOOLS::sqr(s12-ATOOLS::sqr(masses_.hmass))+
       ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
      (ATOOLS::sqr(s12-m_mh2)+m_mh2*m_Gh2);
  }
  for (int n(0);n<2;++n)        GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; sf *= 8./3.; }
  if (j==21) { j=0; sf *= 8./3.; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  epinv_.epinv=epinv2_.epinv2=0.0;
  //if (!m_ishiggs) qqb_ww1jet_v_(p_p,p_msqv);
  //          else 
  gg_hzzg_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  //if (!m_ishiggs) qqb_ww1jet_v_(p_p,p_msqv);
  //           else 
  gg_hzzg_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  //if (!m_ishiggs) qqb_ww1jet_v_(p_p,p_msqv);
  //           else 
  gg_hzzg_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res         * cplfactor * propfactor;
  m_res.IR()     = (res1-res)  * cplfactor * propfactor;
  m_res.IR2()    = (res2-res1) * cplfactor * propfactor;
}

double MCFM_qqb_zz1jet_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_zz1jet_v_Getter,"MCFM_qqb_zz1jet_v")
Virtual_ME2_Base *MCFM_qqb_zz1jet_v_Getter::
operator()(const Process_Info &pi) const
{
  //msg_Out()<<METHOD<<"===================="<<std::endl;
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    //for (int i=0;i<fl.size();i++) msg_Out()<<" "<<fl[i];
    //msg_Out()<<"  ("<<fl.size()<<")."<<std::endl;
    if (fl.size()!=7 && fl.size()!=5) return NULL;
    int pID(0);
    if (fl.size()==5 &&
	(fl[2]==Flavour(kf_Z) && fl[3]==Flavour(kf_Z) && fl[4].Strong())) {
      if (fl[0].IsQuark() && fl[1]==fl[0].Bar()) {
	removebr_.removebr=1;
	pID = 209;
      }
    }
    else if ((fl.size()==7) && fl[6].Strong() &&
	     (fl[2]==fl[3].Bar() && fl[4]==fl[5].Bar() &&
	      fl[2].IsLepton() && fl[2].IsDowntype() &&
	      fl[4].IsLepton() && fl[4].IsDowntype())) {
      if (pi.m_fi.m_ps.size()==2 && fl[0].Strong() && fl[1].Strong()) {
	ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
	ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
	//msg_Out()<<"  check this: "<<fl1<<" & "<<fl2<<std::endl;
	if (fl1==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	    fl2.Strong() && pi.m_fi.m_ps[0].m_ps.size()==2) {
	  ATOOLS::Flavour fl11(pi.m_fi.m_ps[0].m_ps[0].m_fl[0]);
	  ATOOLS::Flavour fl12(pi.m_fi.m_ps[0].m_ps[1].m_fl[0]);
	  if ((fl11==Flavour(kf_Z) && fl12==Flavour(kf_Z))) {
	    removebr_.removebr=0;
	    pID = 209;
	  }
	}
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
      return new MCFM_qqb_zz1jet_v(pi,fl,pID==209);
    }
  }
  return NULL;
}
