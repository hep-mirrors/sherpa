#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;

extern "C" { void qqb_z_v_(double *p,double *msqv); }

namespace MCFM {

  class qqb_z_v: public PHASIC::Virtual_ME2_Base {
  private:

    double *p_p, *p_msqv;
    MODEL::Running_AlphaS *p_as;

  public:

    qqb_z_v(const PHASIC::Process_Info& pi,
	       const ATOOLS::Flavour_Vector& flavs):
      Virtual_ME2_Base(pi,flavs),
      p_as((Running_AlphaS*)s_model->GetScalarFunction("alpha_S"))
    {
      rpa->gen.AddCitation
	(1,"NLO matrix elements from MCFM \\cite{}.");
      p_p = new double[4*MCFM_NMX];
      p_msqv = new double[sqr(2*MCFM_NF+1)];
      m_mode=1;
    }

    ~qqb_z_v()
    {
      delete [] p_p;
      delete [] p_msqv;
    }

    double CallMCFM(const int &i,const int &j)
    {
      qqb_z_v_(p_p,p_msqv);
      return p_msqv[mr(i,j)];
    }
  
    void Calc(const Vec4D_Vector &p)
    {
      for (int n(0);n<2;++n) GetMom(p_p,n,-p[n]);
      for (size_t n(2);n<p.size();++n) GetMom(p_p,n,p[n]);
      long int i(MCFMId(m_flavs[0])), j(MCFMId(m_flavs[1]));
      SetMuR2(m_mur2);
      SetAlphaS((*p_as)(m_mur2));
      epinv_.epinv=epinv2_.epinv2=0.0;
      double res(CallMCFM(i,j));
      epinv_.epinv=1.0;
      double res1(CallMCFM(i,j));
      epinv2_.epinv2=1.0;
      double res2(CallMCFM(i,j));
      const double cfac(4.*9.);
      m_res.Finite()=res*cfac;
      m_res.IR()=(res1-res)*cfac;
      m_res.IR2()=(res2-res1)*cfac;
      m_born=m_res.IR2()/(-2.*qcdcouple_.ason2pi*4./3.*cfac);
    }

    double Eps_Scheme_Factor(const Vec4D_Vector &p)
    {
      return 4.*M_PI;// MSbar scheme
    }

  };// end of class qqb_z_v

}// end of namespace MCFM

using namespace MCFM;

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(qqb_z_v,"qqb_z_v")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,qqb_z_v>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="MCFM") return NULL;
  if (MODEL::s_model->Name()!="SM") return NULL;
  if (!(pi.m_fi.m_nlotype&nlo_type::loop)) return NULL;
  if (pi.m_fi.m_nlocpl[1]!=0.) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size()!=4) return NULL;
  if (fl[0].IsQuark() && fl[1]==fl[0].Bar() &&
      fl[2].IsLepton() && fl[3]==fl[2].Bar()) {
    nproc_.nproc=31;
    chooser_();
    return new qqb_z_v(pi,fl);
  }
  return NULL;
}
