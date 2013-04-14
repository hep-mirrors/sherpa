#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"

namespace MCFM {

  class MCFM_gg_yy: public PHASIC::Tree_ME2_Base {
  private:
    double  m_norm;
    double *p_p, *p_msqv;
    size_t  m_oqcd, m_oew;
    double  m_asfac, m_afac;
    double CouplingFactor(const int oqcd,const int oew) const;
  public:
    MCFM_gg_yy(const PHASIC::Process_Info &pi,
	       const ATOOLS::Flavour_Vector &flavs);
    ~MCFM_gg_yy();
    double Calc(const ATOOLS::Vec4D_Vector &momenta);
    void SetCouplings(const MODEL::Coupling_Map &cpls);
  };

}// end of namespace MCFM

extern "C" { void qqb_gamgam_(double *p,double *msqv); }

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace MCFM;
using namespace MODEL;
using namespace PHASIC;
using namespace ATOOLS;

MCFM_gg_yy::MCFM_gg_yy(const PHASIC::Process_Info& pi,
		       const Flavour_Vector& flavs):
  Tree_ME2_Base(pi,flavs), m_norm(2.0*sqr(2.0*8.0))
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[sqr(2*MCFM_NF+1)];
}

MCFM_gg_yy::~MCFM_gg_yy()
{
  delete [] p_p;
  delete [] p_msqv;
}

void MCFM_gg_yy::SetCouplings(const MODEL::Coupling_Map& cpls)
{
  Tree_ME2_Base::SetCouplings(cpls);
  if (p_aqcd) m_asfac=p_aqcd->Default()/s_model->ScalarFunction("alpha_S");
  if (p_aqed) m_afac=p_aqed->Default()/s_model->ScalarFunction("alpha_QED");
}

double MCFM_gg_yy::CouplingFactor(const int oqcd,const int oew) const
{
  double fac(1.0);
  if (p_aqcd && oqcd) fac*=pow(m_asfac*p_aqcd->Factor(),oqcd);
  if (p_aqed && oew) fac*=pow(m_afac*p_aqed->Factor(),oew);
  return fac;
}

double MCFM_gg_yy::Calc(const Vec4D_Vector &p)
{
  for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
  for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
  scale_.scale=rpa->gen.CplScale();
  scale_.musq=sqr(scale_.scale);
  epinv_.epinv=epinv2_.epinv2=0.0;
  qcdcouple_.as=s_model->ScalarFunction("alpha_S");
  qcdcouple_.gsq=4.*M_PI*qcdcouple_.as;
  qcdcouple_.ason2pi=qcdcouple_.as/(2.*M_PI);
  qcdcouple_.ason4pi=qcdcouple_.as/(4.*M_PI);
  qqb_gamgam_(p_p,p_msqv);
  double res(p_msqv[mr(0,0)]*m_norm);
  return res*CouplingFactor(2,2);
}

extern "C" { void chooser_(); }

DECLARE_TREEME2_GETTER(MCFM_gg_yy,"MCFM_gg_yy")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,MCFM_gg_yy>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype==nlo_type::lo ||
      pi.m_fi.m_nloqcdtype==nlo_type::born) {
    Flavour_Vector fl(pi.ExtractFlavours());
    if (fl.size()!=4) return NULL;
    if (pi.m_fi.m_ps.size()==2) {
      if (fl[0].IsGluon() && fl[1].IsGluon() &&
	  fl[2].IsPhoton() && fl[3].IsPhoton()) {
	if (nproc_.nproc>=0 && nproc_.nproc!=285)
	  THROW(not_implemented,"Only one process allowed in MCFM");
	nproc_.nproc=285;
	chooser_();
	return new MCFM_gg_yy(pi,fl);
      }
    }
  }
  return NULL;
}
