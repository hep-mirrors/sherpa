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

  class MCFM_gg_h_v: public PHASIC::Virtual_ME2_Base {
  private:
    bool    m_plusjet;
    double *p_p, *p_msqv;
    double  m_mh2,m_Gh2,m_mtau2,m_asmh,m_ewcorr;
  public:
    MCFM_gg_h_v(const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs,
		const bool & plusjet);
    ~MCFM_gg_h_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void gg_h_v_(double *p,double *msqv); 
  void gg_hg_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_gg_h_v::MCFM_gg_h_v(const Process_Info& pi,
			 const Flavour_Vector& flavs,
			 const bool & plusjet):
  Virtual_ME2_Base(pi,flavs), m_plusjet(plusjet),
  m_mh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Mass())),
  m_Gh2(ATOOLS::sqr(ATOOLS::Flavour(kf_h0).Width())),
  m_mtau2(ATOOLS::sqr(ATOOLS::Flavour(kf_tau).Mass())),
  m_asmh(MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mh2)),
  m_ewcorr(ATOOLS::sqr(ATOOLS::Flavour(kf_tau).Yuk())/masses_.mbsq/3. *
	   ATOOLS::sqr(1./MODEL::s_model->ScalarConstant(std::string("vev")))/
	   (ewcouple_.gwsq/(4.*ATOOLS::sqr(masses_.wmass))) *
	   ewcouple_.vevsq/
	   ATOOLS::sqr(MODEL::s_model->ScalarConstant(std::string("vev"))))
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_gg_h_v::~MCFM_gg_h_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_gg_h_v::Calc(const Vec4D_Vector &p)
{
  double sf(4.0*9.0/qcdcouple_.ason2pi);
  double cplfactor  = m_ewcorr * ATOOLS::sqr(m_asmh/qcdcouple_.as);
  if (m_plusjet) 
    cplfactor *= 
      MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mur2)/
      qcdcouple_.as;
  double s23((p[2]+p[3]).Abs2());
  double propfactor = 
    (ATOOLS::sqr(s23-ATOOLS::sqr(masses_.hmass))+
     ATOOLS::sqr(masses_.hmass*masses_.hwidth))/
    (ATOOLS::sqr(s23-m_mh2)+m_mh2*m_Gh2);
  for (int n(0);n<2;++n)        GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; sf *= 8./3.; }
  if (j==21) { j=0; sf *= 8./3.; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);

  epinv_.epinv=epinv2_.epinv2=0.0;
  if (!m_plusjet) gg_h_v_(p_p,p_msqv);
              else gg_hg_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  if (!m_plusjet) gg_h_v_(p_p,p_msqv);
             else gg_hg_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  if (!m_plusjet) gg_h_v_(p_p,p_msqv);
             else gg_hg_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res         * cplfactor * propfactor;
  m_res.IR()     = (res1-res)  * cplfactor * propfactor;
  m_res.IR2()    = (res2-res1) * cplfactor * propfactor;
}

double MCFM_gg_h_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_gg_h_v_Getter,"MCFM_gg_h_v")
Virtual_ME2_Base *MCFM_gg_h_v_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if ((fl.size()!=4 && fl.size()!=5) ||
	fl[2]!=fl[3].Bar() ||
	(fl[2]!=ATOOLS::Flavour(kf_tau) && 
	 fl[2]!=ATOOLS::Flavour(kf_tau).Bar())) return NULL;
    int pID(0);
    if (fl.size()==4) {
      if (pi.m_fi.m_ps.size()==1 &&
	  (fl[0].IsGluon() && fl[1].IsGluon())) {
	ATOOLS::Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
	if (flh==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	    ATOOLS::Flavour(kf_tau).Yuk()>0.) pID = 112;
      }
    }
    if (fl.size()==5) {
      if (pi.m_fi.m_ps.size()==2 &&
	  (fl[0].Strong() && fl[1].Strong())) {
	ATOOLS::Flavour flh(pi.m_fi.m_ps[0].m_fl[0]);
	ATOOLS::Flavour flj(pi.m_fi.m_ps[1].m_fl[0]);
	if (flh==ATOOLS::Flavour(kf_h0) && Flavour(kf_h0).IsOn() &&
	    ATOOLS::Flavour(kf_tau).Yuk()>0. &&
	    flj.Strong()) pID = 204;
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
      return new MCFM_gg_h_v(pi,fl,pID==204);
    }
  }
  return NULL;
}
