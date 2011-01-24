#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_qqb_w1jet_v: public PHASIC::Virtual_ME2_Base {
  private:
    double *p_p, *p_msqv;
  public:
    MCFM_qqb_w1jet_v(const PHASIC::Process_Info& pi,
		     const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_w1jet_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { void qqb_w1jet_v_(double *p,double *msqv); }

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_w1jet_v::MCFM_qqb_w1jet_v(const Process_Info& pi,
				   const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs)
{
  p_p    = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_w1jet_v::~MCFM_qqb_w1jet_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_w1jet_v::Calc(const Vec4D_Vector &p)
{
  msg_Out()<<"In "<<METHOD<<"(mu_R^2 = "<<m_mur2<<")."<<std::endl;
  double sf(4.0*9.0);
  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (int n(2);n<5;++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; sf *= 8./3.; }
  if (j==21) { j=0; sf *= 8./3.; }
  scale_.musq    = m_mur2;
  scale_.scale   = sqrt(scale_.musq);
  epinv_.epinv   = epinv2_.epinv2=0.0;
  qqb_w1jet_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv   = 1.0;
  qqb_w1jet_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2 = 1.0;
  qqb_w1jet_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res/qcdcouple_.ason2pi;
  m_res.IR()     = (res1-res)/qcdcouple_.ason2pi;
  m_res.IR2()    = (res2-res1)/qcdcouple_.ason2pi;
  msg_Out()<<"   Results["<<i<<", "<<j<<"]: "<<res2<<" "<<res1<<" "<<res
	   <<" for "<<qcdcouple_.ason2pi<<"."<<std::endl;
}

double MCFM_qqb_w1jet_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_w1jet_v_Getter,"MCFM_qqb_w1jet_v")
Virtual_ME2_Base *MCFM_qqb_w1jet_v_Getter::
operator()(const Process_Info &pi) const
{
  std::cout<<"In "<<METHOD<<":"<<std::endl;
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    msg_Out()<<"   Flavours = "<<fl.size()<<":";
    for (size_t i=0;i<fl.size();i++) msg_Out()<<" "<<fl[i];
    msg_Out()<<"."<<std::endl;
    if (fl.size()!=5) return NULL;
    if (fl[0].Strong() && fl[1].Strong() && fl[4].Strong()) {
      if ((fl[3]==Flavour(kf_nue).Bar() && fl[2]==Flavour(kf_e)) || 
	  (fl[3]==Flavour(kf_numu).Bar() && fl[2]==Flavour(kf_mu)) || 
	  (fl[3]==Flavour(kf_nutau).Bar() && fl[2]==Flavour(kf_tau))) { 
	msg_Out()<<"   --> allowed combination for W-, proc = 16."<<std::endl;
	if (nproc_.nproc>=0) {
	  if (nproc_.nproc!=16)
	    THROW(not_implemented,
		  "Only one process class allowed when using MCFM");
	}
	else {
	  msg_Info()<<"\n";
	  nproc_.nproc=16;
	  chooser_();
	}
      }
      else if ((fl[3]==Flavour(kf_e).Bar() && fl[2]==Flavour(kf_nue)) || 
	       (fl[3]==Flavour(kf_mu).Bar() && fl[2]==Flavour(kf_numu)) || 
	       (fl[3]==Flavour(kf_tau).Bar() && fl[2]==Flavour(kf_nutau))) { 
	msg_Out()<<"   --> allowed combination for W-, proc = 11."<<std::endl;
	if (nproc_.nproc>=0) {
	  if (nproc_.nproc!=11)
	    THROW(not_implemented,
		  "Only one process class allowed when using MCFM");
	}
	else {
	  msg_Info()<<"\n";
	  nproc_.nproc=11;
	  chooser_();
	}
      }
      else {
	THROW(not_implemented,
	      "Check ordering of FS fermions, MCFM wants f\bar f");	
      }
      msg_Out()<<"Will return new MCFM_qqb_w1jet_v."<<std::endl;
      return new MCFM_qqb_w1jet_v(pi,fl);
    }
  }
  return NULL;
}
