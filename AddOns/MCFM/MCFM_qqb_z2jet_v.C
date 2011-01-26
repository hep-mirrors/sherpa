#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_qqb_z2jet_v: public PHASIC::Virtual_ME2_Base {
  private:
    bool    m_neutrinos;
    double *p_p, *p_msqv;
  public:
    MCFM_qqb_z2jet_v(const PHASIC::Process_Info& pi,
		     const ATOOLS::Flavour_Vector& flavs,
		     const bool & neutrinos);
    ~MCFM_qqb_z2jet_v();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { void qqb_z2jet_v_(double *p,double *msqv); }

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_z2jet_v::MCFM_qqb_z2jet_v(const Process_Info& pi,
				   const Flavour_Vector& flavs,
				   const bool & neutrinos):
  Virtual_ME2_Base(pi,flavs), m_neutrinos(neutrinos)
{
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_z2jet_v::~MCFM_qqb_z2jet_v()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_qqb_z2jet_v::Calc(const Vec4D_Vector &p)
{
  double sf(4.0*9.0*(m_neutrinos?1./3.:1.)/qcdcouple_.ason2pi);
  double ason2pi(MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mur2)/
		 (2.*M_PI));
  double asfactor(ATOOLS::sqr(ason2pi/qcdcouple_.ason2pi));
  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (int n(2);n<6;++n) GetMom(p_p,n,p[n]);
  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; sf *= 8./3.; }
  if (j==21) { j=0; sf *= 8./3.; }
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);
  epinv_.epinv=epinv2_.epinv2=0.0;
  qqb_z2jet_v_(p_p,p_msqv);
  double res(p_msqv[mr(i,j)]*sf);
  epinv_.epinv=1.0;
  qqb_z2jet_v_(p_p,p_msqv);
  double res1(p_msqv[mr(i,j)]*sf);
  epinv2_.epinv2=1.0;
  qqb_z2jet_v_(p_p,p_msqv);
  double res2(p_msqv[mr(i,j)]*sf);
  m_res.Finite() = res         * asfactor;
  m_res.IR()     = (res1-res)  * asfactor;
  m_res.IR2()    = (res2-res1) * asfactor;
}

double MCFM_qqb_z2jet_v::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_z2jet_v_Getter,"MCFM_qqb_z2jet_v")
Virtual_ME2_Base *MCFM_qqb_z2jet_v_Getter::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=6) return NULL;
    int pID(0);
    if (fl[0].Strong() && fl[1].Strong() && 
	fl[4].Strong() && fl[5].Strong() &&
	fl[3]==fl[2].Bar()) {
      if (fl[2]==Flavour(kf_e) || fl[2]==Flavour(kf_mu) || fl[2]==Flavour(kf_tau)) {
	pID = 44;
      }
      else if (fl[2]==Flavour(kf_nue) || fl[2]==Flavour(kf_numu) || fl[2]==Flavour(kf_nutau)) {
	pID = 46;
      }
    } 
    if (pID>0) {
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      else {
	nproc_.nproc=pID;
	chooser_();
      }
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_qqb_z2jet_v(pi,fl,pID==46);
    }
  }
  return NULL;
}
