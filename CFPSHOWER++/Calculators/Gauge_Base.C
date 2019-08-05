#include "CFPSHOWER++/Calculators/Gauge_Base.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::Gauge_Base
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"


using namespace CFPSHOWER;
using namespace ATOOLS;

Gauge_Base::Gauge_Base(const Kernel_Info & info) :
  m_type(info.Type()),
  m_tagsequence(info.TagSequence()),
  p_alphaS(info.GetAlphaS()), p_alpha(info.GetAlpha()),
  m_CF(4./3.), m_CA(3.), m_TR(1./2.), m_zeta3(1.202056903159594),
  m_charge(1.), 
  m_muR2factor(info.MuR2Factor()), m_asfactor(info.AsFactor()), 
  m_orderA(1), m_orderB(1),
  m_name("undefined")
{
  // KFactor: 1 = include A2
  // KFactor: 2 = include B2
  // KFactor: 4 = include A3
  switch (info.KFactor()) {
  case 7: m_orderA = 3; m_orderB = 2; break;
  case 5: m_orderA = 3; m_orderB = 1; break;
  case 3: m_orderA = 2; m_orderB = 2; break;
  case 1: m_orderA = 2; m_orderB = 1; break;
  case 0: m_orderA = 1; m_orderB = 1; break;
  default:
    msg_Error()<<"Error in "<<METHOD
	       <<": invalid higher order setting, will use LO.\n";
    break;
  }
  if (p_alphaS) {
    m_K1max  = K1(3.);
    m_K2max  = K2(3.);
  }
}

double Gauge_Base::operator()(const Splitting & split) {
  return (*p_alphaS)(m_muR2factor*Scale(split));
}

const double Gauge_Base::OverEstimate(const Splitting & split) const {
  if (p_alphaS) return (*p_alphaS)(m_muR2factor*split.t0());
}

const double Gauge_Base::Scale(const Splitting & split) const {
  return split.t();
}

const double Gauge_Base::Beta0(const double & NF) const {
  return 11./6.*m_CA - 2./3.*m_TR*NF;
}

const double Gauge_Base::Beta1(const double & NF) const {
  return 17./6.*sqr(m_CA) - (5./3.*m_CA+m_CF)*m_TR*NF;
}

const double Gauge_Base::NF(const double & q2) const {
  return p_alphaS->Nf(q2);
}
  
const double Gauge_Base::K(const Splitting & split) const {
  return K(Scale(split));
}

const double Gauge_Base::KMax(const Splitting & split) const {
  if (m_orderA<2) return 0.;
  double alphaS = (*p_alphaS)(split.t0())/(2.*M_PI);
  if (m_orderA<3) return alphaS * m_K1max;
  return alphaS * m_K1max + sqr(alphaS) * m_K2max;
}

const double Gauge_Base::K(const double & q2) const {
  // Coefficients K1 & K2 agree with DIRE.
  if (m_orderA<2) return 0.;
  double nf = NF(q2), alphaS = (*p_alphaS)(q2)/(2.*M_PI);
  if (m_orderA<3) return alphaS * K1(nf);
  return alphaS * K1(nf) + sqr(alphaS) * K2(nf);
}

const double Gauge_Base::K1(const double & NF) const {
  return m_CA*(67./18.-sqr(M_PI)/6.)-10./9.*m_TR*NF;
}

const double Gauge_Base::K2(const double & NF) const {
  return (sqr(m_CA)     * (245./6.-134./27.*sqr(M_PI)+11./45.*pow(M_PI,4) +
			   22./3.*m_zeta3)
	  +m_CA*m_TR*NF * (-418./27.+40./27.*sqr(M_PI)-56./3.*m_zeta3)
	  +m_CF*m_TR*NF * (-55./3.+16.*m_zeta3)-16./27.*sqr(m_TR*NF))/4.;
}

