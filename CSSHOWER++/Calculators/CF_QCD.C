#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Exception.H"

namespace CSSHOWER {
  
  const double s_Nc = 3.;
  const double s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  const double s_CA = s_Nc;
  const double s_TR = 1./2.;

  class CF_QCD: public SF_Coupling {
  protected:

    MODEL::Running_AlphaS *p_cpl;

    double m_q;

  public:

    inline CF_QCD(const SF_Key &key):
      SF_Coupling(key) 
    {
      if (key.p_v->in[0].StrongCharge()==8 &&
	  key.p_v->in[1].StrongCharge()==8 &&
	  key.p_v->in[2].StrongCharge()==8) m_q=s_CA;
      else m_q=(key.p_v->in[0].StrongCharge()==8)?s_TR:s_CF;
    }

    bool SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int mode);
    bool AllowSpec(const ATOOLS::Flavour &fl);

    double CplFac(const double &scale) const;

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

bool CF_QCD::SetCoupling(MODEL::Model_Base *md,const double &k0sq,
			 const double &isfac,const double &fsfac)
{
  p_cpl=(MODEL::Running_AlphaS*)md->GetScalarFunction("alpha_S");
  m_cplfac=((m_type/10==1)?fsfac:isfac)/CplFac(k0sq);
  m_cplmax.push_back((*p_cpl)(k0sq)*m_q);
  m_cplmax.push_back(0.0);
  return true;
}

double CF_QCD::Coupling(const double &scale,const int mode)
{
  if (mode!=0) return 0.0;
  double cpl=(*p_cpl)(CplFac(scale)*scale)*m_q;
  if (cpl>m_cplmax.front()) return m_cplmax.front();
  return cpl;
}

bool CF_QCD::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (abs(fl.StrongCharge())==3) {
    switch (m_type) {
    case cstp::FF: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::FI: 
      if (abs(p_lf->FlA().StrongCharge())==3)
	return p_lf->FlA().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::IF: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==fl.StrongCharge();
      break;
    case cstp::II: 
      if (abs(p_lf->FlB().StrongCharge())==3)
	return p_lf->FlB().StrongCharge()==-fl.StrongCharge();
      break;
    case cstp::none: abort();
    }
  }
  return fl.Strong();
}

double CF_QCD::CplFac(const double &scale) const
{
  if (m_kfmode==0) return m_cplfac;
  double nf=p_cpl->Nf(scale);
  double kfac=exp(-(67.0-3.0*sqr(M_PI)-5.0*nf)/(33.0-2.0*nf));
  return m_cplfac*kfac;
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
  if (MODEL::s_model->Name().find("MSSM")==std::string::npos) return NULL;
  std::string sgtag("{"+Flavour(kf_Gluino).IDName()+"}");
  new CF_QCD_Getter(sgtag+sgtag+gtag);
  new CF_QCD_Getter(sgtag+gtag+sgtag);
  new CF_QCD_Getter(gtag+sgtag+sgtag);
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(1000000+i));
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    new CF_QCD_Getter(gtag+qtag+qbtag);
    new CF_QCD_Getter(qbtag+qbtag+gtag);
    new CF_QCD_Getter(qtag+qtag+gtag);
  }
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(2000000+i));
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
