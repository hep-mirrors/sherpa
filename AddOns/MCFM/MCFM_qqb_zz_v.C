#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_qqb_zz_v: public PHASIC::Virtual_ME2_Base {
  private:
    double *p_p, *p_msqv;
  public:
    MCFM_qqb_zz_v(const PHASIC::Process_Info& pi,
		  const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_zz_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { void qqb_zz_v_(double *p,double *msqv); }

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_zz_v::MCFM_qqb_zz_v(const Process_Info& pi,
			   const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs)
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_zz_v::~MCFM_qqb_zz_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_zz_v::Calc(const Vec4D_Vector &p)
{
  const double sf(4.0*9.0/qcdcouple_.ason2pi);
  for (size_t n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);
  epinv_.epinv=epinv2_.epinv2=0.0;
  qqb_zz_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  qqb_zz_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  qqb_zz_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
}

double MCFM_qqb_zz_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_zz_v_Getter,"MCFM_qqb_zz_v")
Virtual_ME2_Base *MCFM_qqb_zz_v_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=6) return NULL;
    if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[3]==fl[2].Bar() && fl[5]==fl[4].Bar() &&
	(fl[2]==Flavour(kf_e) || fl[2]==Flavour(kf_mu) || fl[2]==Flavour(kf_tau)) &&
	(fl[4]==Flavour(kf_e) || fl[4]==Flavour(kf_mu) || fl[4]==Flavour(kf_tau))) {
      double yuk(MODEL::s_model->ScalarConstant("Yukawa_"+ToString(fl[0])));
      if (yuk>0.0 && Flavour(kf_h0).IsOn()) return NULL;
      int pid=86;
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pid)
	  THROW(not_implemented,"Only one process class allowed when using MCFM");
      }
      else {
	msg_Info()<<"\n";
	nproc_.nproc=pid;
	chooser_();
      }
      return new MCFM_qqb_zz_v(pi,fl);
    }
    if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
	fl[4]==fl[2].Bar() && fl[5]==fl[3].Bar() && fl[2]!=fl[3] &&
	(fl[2]==Flavour(kf_e) || fl[2]==Flavour(kf_mu) || fl[2]==Flavour(kf_tau)) &&
	(fl[3]==Flavour(kf_e) || fl[3]==Flavour(kf_mu) || fl[3]==Flavour(kf_tau))) {
      int pid=81;
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pid)
	  THROW(not_implemented,"Only one process class allowed when using MCFM");
      }
      else {
	msg_Info()<<"\n";
	nproc_.nproc=pid;
	chooser_();
      }
      return new MCFM_qqb_zz_v(pi,fl);
    }
  }
  return NULL;
}
