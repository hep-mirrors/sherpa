#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"

namespace CSSHOWER {
  
  class CF_GGH: public SF_Coupling {
  public:

    inline CF_GGH(const SF_Key &key): SF_Coupling(key) {}

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace ATOOLS;

bool CF_GGH::SetCoupling(MODEL::Model_Base *md,
			 const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  double vev(std::abs(md->ComplexConstant("cvev")));
  double asggh(md->ScalarFunction(std::string("alpha_S"),
				  sqr(Flavour(kf_h0).Mass())));
  double cpl(asggh/(2.0*M_PI*vev));
  DEBUG_VAR(cpl);
  m_cplfac=1.0;
  m_cplmax.push_back(cpl*cpl);
  return true;
}

double CF_GGH::Coupling(const double &scale,const int pol)
{
  if (pol>0) return 0.0;
  return m_last = m_cplmax.front();
}

bool CF_GGH::AllowSpec(const ATOOLS::Flavour &fl) 
{
  return true;
}

namespace CSSHOWER {

DECLARE_CPL_GETTER(CF_GGH_Getter);

SF_Coupling *CF_GGH_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_GGH(args);
}

void CF_GGH_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"ggh coupling";
}

}

DECLARE_GETTER(CF_GGH_Getter,"SF_GGH_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_GGH_Getter>::
operator()(const SFC_Filler_Key &key) const
{
  if (!Flavour(kf_h0).IsOn()) return NULL;
  std::string gtag("{"+Flavour(kf_gluon).IDName()+"}");
  std::string htag("{"+Flavour(kf_h0).IDName()+"}");
  key.p_gets->push_back(new CF_GGH_Getter(gtag+gtag+htag));
  key.p_gets->push_back(new CF_GGH_Getter(gtag+htag+gtag));
  key.p_gets->push_back(new CF_GGH_Getter(htag+gtag+gtag));
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_GGH_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"ggh coupling filler";
}
