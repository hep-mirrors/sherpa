#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {
  class MCFM_qqb_QQb: public PHASIC::Virtual_ME2_Base {
  private:
    int     m_pID;
    bool    m_swapped;
    double  m_normcpl, m_normcorr;
    double *p_p, *p_msqv;

    double CallMCFM(const int & i,const int & j);
  public:
    MCFM_qqb_QQb(const int & pID,const bool & swapped, 
		const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs);
    ~MCFM_qqb_QQb();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };

}// end of namespace MCFM

extern "C" { 
  void qqb_qqb_v_(double *p,double *msqv); 
}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_qqb_QQb::MCFM_qqb_QQb(const int & pID,const bool & swapped,
			 const PHASIC::Process_Info& pi,
			 const Flavour_Vector& flavs):
  Virtual_ME2_Base(pi,flavs), m_pID(pID), m_swapped(swapped),
  m_normcpl(4.*9./qcdcouple_.ason2pi)
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
  m_drmode=m_mode=1;
}

MCFM_qqb_QQb::~MCFM_qqb_QQb()
{
  delete [] p_p;
  delete [] p_msqv;
}

double MCFM_qqb_QQb::CallMCFM(const int & i,const int & j) {
  qqb_qqb_v_(p_p,p_msqv);
  return p_msqv[mr(i,j)];
}

void MCFM_qqb_QQb::Calc(const Vec4D_Vector &p)
{
  scale_.musq=m_mur2;
  scale_.scale=sqrt(scale_.musq);
  double corrfactor(m_normcpl);
  double as(MODEL::s_model->ScalarFunction(std::string("alpha_S"),m_mur2));
  if (!m_flavs[0].IsQuark()) {
    corrfactor *= 64./9.;
  }
  corrfactor *= pow(4.*M_PI*as/qcdcouple_.gsq,2);
  
  msg_Debugging()<<METHOD<<" with corr = "<<corrfactor<<" for "
	   <<m_flavs[0]<<" & "<<m_flavs[1]<<".\n";

  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  if (m_pID==157) {
    for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  }
  else {
    for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  }

  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }
  if (j==21) { j=0; }

  epinv_.epinv=epinv2_.epinv2=0.0;
  double res(CallMCFM(i,j)  * corrfactor);
  epinv_.epinv=1.0;
  double res1(CallMCFM(i,j) * corrfactor);
  epinv2_.epinv2=1.0;
  double res2(CallMCFM(i,j) * corrfactor);

  msg_Debugging()<<"   --> "<<res<<" "<<res1<<" "<<res2<<" "
	   <<"for corrfactor/as = "<<corrfactor<<", "<<as<<".\n";
  m_res.Finite() = res;
  m_res.IR()     = (res1-res);
  m_res.IR2()    = (res2-res1);
}

double MCFM_qqb_QQb::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_QQb_Getter,"MCFM_qqb_QQb")
Virtual_ME2_Base *MCFM_qqb_QQb_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    if (fl.size()!=4)                                   return NULL;
    int pID(0);
    bool swapped(false);
    if (pi.m_fi.m_ps.size()==2) {
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
      if (fl1==Flavour(kf_t) && fl2==Flavour(kf_t).Bar()) {
	pID = 157;
        zerowidth_.zerowidth=true;
	swapped=false;
      }
    }
    else {
      return NULL;
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
      return new MCFM_qqb_QQb(pID,swapped,pi,fl);
    }
  }
  return NULL;
}
