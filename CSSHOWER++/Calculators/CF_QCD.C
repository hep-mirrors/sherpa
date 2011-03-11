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
      if (key.m_type==cstp::FF || key.m_type==cstp::FI) {
	if (key.p_v->in[0].StrongCharge()==8) m_q/=2.0;
      }
      else {
	if (key.m_mode==0) {
	  if (key.p_v->in[1].StrongCharge()==8) m_q/=2.0;
	}
	else {
	  if (key.p_v->in[2].StrongCharge()==8) m_q/=2.0;
	}
      }
    }

    bool SetCoupling(MODEL::Model_Base *md,const double &k0sq,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
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
  m_cplfac=((m_type/10==1)?fsfac:isfac);
  m_cplmax.push_back((*p_cpl)(k0sq)*m_q);
  m_cplmax.push_back(0.0);
  return true;
}

double CF_QCD::Coupling(const double &scale,const int pol)
{
  if (pol!=0) return 0.0;
  double scl(CplFac(scale)*scale);
  double cpl=(*p_cpl)(scl)*m_q;
  if (cpl>m_cplmax.front()) return m_cplmax.front();
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"as weight kt = "<<sqrt(CplFac(scale))<<" * "
		 <<sqrt(scale)<<", \\alpha_s("<<sqrt(scl)<<") = "
		 <<(*p_cpl)(scl)<<", m_q = "<<m_q<<"\n";
#endif
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
  double kfac=exp(-(67.0-3.0*sqr(M_PI)-10.0/3.0*nf)/(33.0-2.0*nf));
  return m_cplfac*kfac;
}

namespace CSSHOWER {

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
	       void,SFC_Filler_Key);

void *CF_QCD_Filler::operator()
  (const SFC_Filler_Key &key) const
{
  if (!Flavour(kf_gluon).IsOn()) return NULL;
  std::string gtag("{"+Flavour(kf_gluon).IDName()+"}");
  key.p_gets->push_back(new CF_QCD_Getter(gtag+gtag+gtag));
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)i);
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  if (MODEL::s_model->Name().find("MSSM")==std::string::npos) return NULL;
  std::string sgtag("{"+Flavour(kf_Gluino).IDName()+"}");
  key.p_gets->push_back(new CF_QCD_Getter(sgtag+sgtag+gtag));
  key.p_gets->push_back(new CF_QCD_Getter(sgtag+gtag+sgtag));
  key.p_gets->push_back(new CF_QCD_Getter(gtag+sgtag+sgtag));
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(1000000+i));
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  for (int i(1);i<=6;++i) {
    Flavour f((kf_code)(2000000+i));
    if (!f.IsOn()) continue;
    std::string qtag("{"+f.IDName()+"}");
    std::string qbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_QCD_Getter(gtag+qtag+qbtag));
    key.p_gets->push_back(new CF_QCD_Getter(qbtag+qbtag+gtag));
    key.p_gets->push_back(new CF_QCD_Getter(qtag+qtag+gtag));
  }
  return NULL;
}

void CF_QCD_Filler::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"qcd coupling filler";
}

}
