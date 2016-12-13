#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

namespace CSSHOWER {
  
  class CF_HV: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;

    double m_cplfac, m_q;

  public:

    inline CF_HV(const SF_Key &key):
      SF_Coupling(key) 
    {
      if (key.p_v->in[0].StrongCharge()==8 &&
	  key.p_v->in[1].StrongCharge()==8 &&
	  key.p_v->in[2].StrongCharge()==8) m_q=0;
      else m_q=(key.p_v->in[0].StrongCharge()==8)?1:2;
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
    
    bool   SetCoupling(MODEL::Model_Base *md,
		       const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
    bool   AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

bool CF_HV::SetCoupling(MODEL::Model_Base *md,
			const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  double CF = double(md->ScalarConstant(std::string("CF")));
  double CA = double(md->ScalarConstant(std::string("CA")));
  double TR = double(md->ScalarConstant(std::string("TR")));

  p_cpl=md->GetScalarFunction("alpha_HV");
  switch (int(m_q)) {
  case 0: m_q=CA;break; 
  case 1: m_q=TR;break; 
  case 2: m_q=CF;break;
  }
  m_cplfac=((m_type/10==1)?fsfac:isfac);
  double scale((m_type/10==1)?k0sqf:k0sqi);
  double scl(CplFac(scale)*scale);
  m_cplmax.push_back((*p_cpl)(scl)*m_q);
  std::cout<<" cpl max HV "<<m_cplmax.back()<<" "
	   <<k0sqi<<"/"<<k0sqf<<" "<<m_q <<std::endl; 
  m_cplmax.push_back(0.0);
  return true;
}

double CF_HV::Coupling(const double &scale,const int pol)
{
  if (pol!=0) return 0.0;
  double scl(CplFac(scale)*scale);
  return (*p_cpl)(scl)*m_q;
}

bool CF_HV::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return (fl.Strong()&&fl.Kfcode()>9900000);
}

namespace CSSHOWER {

DECLARE_CPL_GETTER(CF_HV_Getter);

SF_Coupling *CF_HV_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_HV(args);
}

void CF_HV_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"HV gauge coupling";
}

}

DECLARE_GETTER(CF_HV_Getter,"SF_HV_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_HV_Getter>::
operator()(const SFC_Filler_Key &key) const
{
  if (key.p_md->Name()!=std::string("SM+HiddenValley")) return NULL;
  if (!Flavour(9900021).IsOn()) return NULL;
  std::string gDtag("{"+Flavour(9900021).IDName()+"}");
  key.p_gets->push_back(new CF_HV_Getter(gDtag+gDtag+gDtag));
  for (int i(1);i<3;++i) {
    Flavour f((kf_code)9900000+i);
    if (!f.IsOn()) continue;
    std::string qDtag("{"+f.IDName()+"}");
    std::string qDbtag ("{"+f.Bar().IDName()+"}");
    key.p_gets->push_back(new CF_HV_Getter(gDtag+qDtag+qDbtag));
    key.p_gets->push_back(new CF_HV_Getter(qDbtag+qDbtag+gDtag));
    key.p_gets->push_back(new CF_HV_Getter(qDtag+qDtag+gDtag));
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_HV_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"HV coupling filler";
}
