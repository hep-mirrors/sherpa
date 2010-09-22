#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_qqb_ww_v: public PHASIC::Virtual_ME2_Base {
  private:
    double *p_p, *p_msqv;
  public:
    MCFM_qqb_ww_v(const PHASIC::Process_Info& pi,
		 const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_ww_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { void qqb_ww_v_(double *p,double *msqv); }

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_ww_v::MCFM_qqb_ww_v(const Process_Info& pi,
			   const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs)
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
  const double sf(4.0*9.0);
  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (int n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);
  epinv_.epinv=epinv2_.epinv2=0.0;
  qqb_ww_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  qqb_ww_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  qqb_ww_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res/qcdcouple_.ason2pi;
  m_res.IR()     = (res1-res)/qcdcouple_.ason2pi;
  m_res.IR2()    = (res2-res1)/qcdcouple_.ason2pi;
}

double MCFM_qqb_ww_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_ww_v_Getter,"MCFM_qqb_ww_v")
Virtual_ME2_Base *MCFM_qqb_ww_v_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl[0].IsQuark() && fl[1]==fl[0].Bar()) {
      if (fl.size()==4 &&
	  (fl[2]==Flavour(kf_Wplus) && fl[3]==Flavour(kf_Wplus).Bar())) {
	removebr_.removebr=1;
      }
    }
    else if (fl.size()==6 &&
	     ((fl[2]==Flavour(kf_nue)   && fl[3]==Flavour(kf_e).Bar()) || 
	      (fl[2]==Flavour(kf_numu)  && fl[3]==Flavour(kf_mu).Bar()) || 
	      (fl[2]==Flavour(kf_nutau) && fl[3]==Flavour(kf_tau).Bar())) &&
	     ((fl[4]==Flavour(kf_e)     && fl[2]==Flavour(kf_nue).Bar()) || 
	      (fl[4]==Flavour(kf_mu)    && fl[2]==Flavour(kf_numu).Bar()) || 
	      (fl[4]==Flavour(kf_tau)   && fl[2]==Flavour(kf_nutau).Bar()))) {
	removebr_.removebr=0;
    }
    else {
      THROW(not_implemented,"Check ordering of FS fermions, MCFM wants f\bar f");
    }
    int pid=61;
    zerowidth_.zerowidth=true;
    if (nproc_.nproc>=0) {
      if (nproc_.nproc!=pid)
	THROW(not_implemented,"Only one process class allowed when using MCFM");
    }
    msg_Info()<<"\n";
    nproc_.nproc=pid;
    chooser_();
    msg_Info()<<"   Return MCFM interface for pid = "<<pid<<"."<<std::endl;
    PRINT_VAR(qcdcouple_.ason2pi<<" "<<scale_.musq);
    return new MCFM_qqb_ww_v(pi,fl);
  }
  return NULL;
}
