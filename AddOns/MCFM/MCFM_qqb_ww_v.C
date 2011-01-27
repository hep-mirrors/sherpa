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



  class MCFM_qqb_ww_v: public PHASIC::Virtual_ME2_Base {
  private:
    bool    m_ishiggs;
    double *p_p, *p_msqv;
    double  m_mh2,m_Gh2,m_asmh,m_ewcorr;
  public:
    MCFM_qqb_ww_v(const PHASIC::Process_Info& pi,
		  const ATOOLS::Flavour_Vector& flavs,
		  const bool & ishiggs);
    ~MCFM_qqb_ww_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_ww_v_(double *p,double *msqv); 
  void qqb_hww_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_ww_v::MCFM_qqb_ww_v(const Process_Info& pi,
			     const Flavour_Vector& flavs,
			     const bool & ishiggs):
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

MCFM_qqb_ww_v::~MCFM_qqb_ww_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_ww_v::Calc(const Vec4D_Vector &p)
{
  double sf(4.0*(m_ishiggs?64.0:9.0)/qcdcouple_.ason2pi);
  double cplfactor(1.),propfactor(1.);
  if (m_ishiggs) {
    cplfactor   = m_ewcorr * ATOOLS::sqr(m_asmh/qcdcouple_.as);
    double s12((p[0]+p[1]).Abs2());
    propfactor = 
      (ATOOLS::sqr(s12-ATOOLS::sqr(masses_.hmass))+
       ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
      (ATOOLS::sqr(s12-m_mh2)+m_mh2*m_Gh2);
  }
  for (int n(0);n<2;++n)        GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  epinv_.epinv=epinv2_.epinv2=0.0;
  if (!m_ishiggs) qqb_ww_v_(p_p,p_msqv);
             else qqb_hww_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  if (!m_ishiggs) qqb_ww_v_(p_p,p_msqv);
             else qqb_hww_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  if (!m_ishiggs) qqb_ww_v_(p_p,p_msqv);
             else qqb_hww_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res         * cplfactor * propfactor;
  m_res.IR()     = (res1-res)  * cplfactor * propfactor;
  m_res.IR2()    = (res2-res1) * cplfactor * propfactor;
}

double MCFM_qqb_ww_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_ww_v_Getter,"MCFM_qqb_ww_v")
Virtual_ME2_Base *MCFM_qqb_ww_v_Getter::operator()(const Process_Info &pi) const
{
  //  msg_Out()<<METHOD<<"===================="<<std::endl;
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=6 && fl.size()!=4) return NULL;
    int pID(0);
    if (fl.size()==4 &&
	(fl[2]==Flavour(kf_Wplus) && fl[3]==Flavour(kf_Wplus).Bar())) {
      if (fl[0].IsQuark() && fl[1]==fl[0].Bar()) {
	removebr_.removebr=1;
	pID = 61;
      }
    }
    else if ((fl.size()==6) &&
	     ((((fl[2]==Flavour(kf_nue)   && fl[3]==Flavour(kf_e).Bar()) || 
		(fl[3]==Flavour(kf_nue)   && fl[2]==Flavour(kf_e).Bar()) || 
		(fl[2]==Flavour(kf_numu)  && fl[3]==Flavour(kf_mu).Bar()) || 
		(fl[3]==Flavour(kf_numu)  && fl[2]==Flavour(kf_mu).Bar()) || 
		(fl[2]==Flavour(kf_nutau) && fl[3]==Flavour(kf_tau).Bar()) ||
		(fl[3]==Flavour(kf_nutau) && fl[2]==Flavour(kf_tau).Bar())) &&
	       ((fl[4]==Flavour(kf_e)     && fl[5]==Flavour(kf_nue).Bar()) || 
		(fl[5]==Flavour(kf_e)     && fl[4]==Flavour(kf_nue).Bar()) || 
		(fl[4]==Flavour(kf_mu)    && fl[5]==Flavour(kf_numu).Bar()) || 
		(fl[5]==Flavour(kf_mu)    && fl[4]==Flavour(kf_numu).Bar()) || 
		(fl[4]==Flavour(kf_tau)   && fl[5]==Flavour(kf_nutau).Bar()) ||
		(fl[5]==Flavour(kf_tau)   && fl[4]==Flavour(kf_nutau).Bar()))) ||
	      (((fl[4]==Flavour(kf_nue)   && fl[5]==Flavour(kf_e).Bar()) || 
		(fl[5]==Flavour(kf_nue)   && fl[4]==Flavour(kf_e).Bar()) || 
		(fl[4]==Flavour(kf_numu)  && fl[5]==Flavour(kf_mu).Bar()) || 
		(fl[5]==Flavour(kf_numu)  && fl[4]==Flavour(kf_mu).Bar()) || 
		(fl[4]==Flavour(kf_nutau) && fl[5]==Flavour(kf_tau).Bar()) ||
		(fl[5]==Flavour(kf_nutau) && fl[4]==Flavour(kf_tau).Bar())) &&
	       ((fl[2]==Flavour(kf_e)     && fl[3]==Flavour(kf_nue).Bar()) || 
		(fl[3]==Flavour(kf_e)     && fl[2]==Flavour(kf_nue).Bar()) || 
		(fl[2]==Flavour(kf_mu)    && fl[3]==Flavour(kf_numu).Bar()) || 
		(fl[3]==Flavour(kf_mu)    && fl[2]==Flavour(kf_numu).Bar()) || 
		(fl[2]==Flavour(kf_tau)   && fl[3]==Flavour(kf_nutau).Bar()) ||
		(fl[3]==Flavour(kf_tau)   && fl[2]==Flavour(kf_nutau).Bar()))))) {
      if (pi.m_fi.m_ps.size()==2 &&
	  (fl[0].IsQuark() && fl[1]==fl[0].Bar())) {
	ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]),fl2(pi.m_fi.m_ps[1].m_fl[0]);
	if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	    (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	  removebr_.removebr=0;
	  pID = 61;
	}
      }
      else if (pi.m_fi.m_ps.size()==1 && (fl[0].IsGluon() && fl[1]==fl[0])) {
	ATOOLS::Flavour fl0(pi.m_fi.m_ps[0].m_fl[0]);
	if (fl0==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	    pi.m_fi.m_ps[0].m_ps.size()==2) {
	  ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_ps[0].m_fl[0]);
	  ATOOLS::Flavour fl2(pi.m_fi.m_ps[0].m_ps[1].m_fl[0]);
	  if ((fl1==Flavour(kf_Wplus) && fl2==Flavour(kf_Wplus).Bar()) ||
	      (fl2==Flavour(kf_Wplus) && fl1==Flavour(kf_Wplus).Bar())) {
	    removebr_.removebr=0;
	    pID = 113;
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
      return new MCFM_qqb_ww_v(pi,fl,pID==113);
    }
  }
  return NULL;
}
