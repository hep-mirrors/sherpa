#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"

namespace CSSHOWER {
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;

    double m_cplfac, m_q;

  public:

    inline CF_QCD(const SF_Key &key):
      SF_Coupling(key) 
    {
      if (key.p_v->in[0].StrongCharge()==8 &&
	  key.p_v->in[1].StrongCharge()==8 &&
	  key.p_v->in[2].StrongCharge()==8) m_q=s_CA;
      else m_q=(key.p_v->in[0].StrongCharge()==8)?s_TR:s_CF;
    }

    void SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int mode);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

void CF_QCD::SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			 const double &isfac,const double &fsfac)
{
  m_cplfac=(m_type/10==1)?fsfac:isfac;
  p_cpl=md->GetScalarFunction("alpha_S");
  m_cplmax.push_back((*p_cpl)(m_cplfac*k0sq)*m_q);
  m_cplmax.push_back(0.0);
}

double CF_QCD::Coupling(const double &scale,const int mode)
{
  if (mode!=0) return 0.0;
  return (*p_cpl)(m_cplfac*scale)*m_q;
}

bool CF_QCD::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return fl.Strong();
}

DECLARE_CPL_GETTER(CF_QCD_Getter);

SF_Coupling *CF_QCD_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_QCD(args);
}

void CF_QCD_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"strong coupling";
}

DECLARE_GETTER(CF_QCD_Filler,"SF_QCD_Fill",
	       void,const MODEL::Model_Base *);

void *CF_QCD_Filler::operator()
  (const MODEL::Model_Base *const &model) const
{
  if (!Flavour(kf_gluon).IsOn()) return NULL;
  std::string gtag("{"+Flavour(kf_gluon).IDName()+"}");
  new CF_QCD_Getter(gtag+gtag+gtag);
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)i);
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    new CF_QCD_Getter(gtag+qtag+qbtag);
    new CF_QCD_Getter(qbtag+qbtag+gtag);
    new CF_QCD_Getter(qtag+qtag+gtag);
  }
  return NULL;
}

void CF_QCD_Filler::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd coupling filler";
}
