#include "HADRONS++/PS_Library/Resonance.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


using namespace HADRONS;
using namespace ATOOLS;
using namespace std;


Resonance_Base::Resonance_Base(const Resonance_Parameters & params) :
  m_inflav(params.m_inflav), m_type(params.m_type),
  m_OSmass(params.m_OSmass<0.   ? m_inflav.Mass(true) : params.m_OSmass), m_OSmass2(m_OSmass*m_OSmass),
  m_OSwidth(params.m_OSwidth<0. ? m_inflav.Width() : params.m_OSwidth),   m_OSwidth2(m_OSwidth*m_OSwidth),
  m_phase(params.m_phase), m_threshold(0.), m_threshold2(0.), m_lambda(0.),
  m_exponent(m_inflav.IsVector() ? 3 : (m_inflav.IsScalar() ? 1 : 0)) 
{
  m_name = std::string("R_")+m_inflav.IDName();
  if (m_type==resonance_type::bespoke)      m_name += std::string("bespoke");
  else if (m_type==resonance_type::running) m_name += std::string("running");
  else                                      m_name += std::string("fixed");
  m_decmasses2.resize(params.m_outflavs.size(),0.);
  for (size_t i=0;i<params.m_outflavs.size();i++) {
    double mass     = params.m_outflavs[i].Mass(true);
    m_threshold    += mass;
    m_decmasses2[i] = ATOOLS::sqr(mass);
  }
  m_threshold2  = ATOOLS::sqr(m_threshold);
  m_lambda      = Lambda(m_OSmass2,m_decmasses2[0],m_decmasses2[1]);
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

RunningWidth_Resonance::RunningWidth_Resonance(const Resonance_Parameters & params) :
  Resonance_Base(params)
{}

double RunningWidth_Resonance::CalculateWidth(const double & s) {
  if (s>m_threshold2) {
    if (dabs(1.-s/m_OSmass2)<1.e-12) return m_OSwidth;
    if (m_decmasses2.size()==2) {
      return ( m_OSwidth *
	       pow(Lambda(s,m_decmasses2[0],m_decmasses2[1])/(2.*sqrt(s)) /
		   (m_lambda / (2.*m_OSmass)), m_exponent)
	       );
    }
  }
  return 0.;
}

Complex RunningWidth_Resonance::BreitWigner(const double & s) {
  double MW = sqrt(s)*CalculateWidth(s);
  return m_OSmass2/Complex( m_OSmass2-s, -MW );
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////


GS_Resonance::GS_Resonance(const Resonance_Parameters & params) : 
  Resonance_Base(params),
  m_m_pi(Flavour(kf_pi_plus).Mass(true))
{
  double mpi2 = sqr(m_m_pi);
  m_d   = ( 3/M_PI * mpi2/sqr(m_lambda) * log((m_OSmass+2.*m_lambda)/(2.*m_m_pi)) +
	    m_OSmass/(2.*M_PI*m_lambda) * (1. - 2.*mpi2/sqr(m_lambda)) );
  m_h0  = h(m_OSmass2);
  m_dh0 = m_h0 * (1./(8.*sqr(m_lambda))-1./(2.*m_OSmass2));
}

double GS_Resonance::CalculateWidth(const double & s) {
  if (s>m_threshold2) {
    if (dabs(1.-s/m_OSmass2)<1.e-12) return m_OSwidth;
    if (m_decmasses2.size()==2) {
      return ( m_OSwidth *
	       pow(Lambda(s,m_decmasses2[0],m_decmasses2[1])/sqrt(s) /
		   (m_lambda / m_OSmass), m_exponent)
	       );
    }
  }
  return m_OSwidth;
}

Complex GS_Resonance::BreitWigner(const double & s) {
  double Gamma = CalculateWidth(s);
  return ( (m_OSmass2+m_d*m_OSmass*m_OSwidth)/
	   Complex(m_OSmass2-s+f(s), -sqrt(s)*Gamma) );
}

double GS_Resonance::h(const double & s) {
  double ks = Lambda(s,m_decmasses2[0],m_decmasses2[1]), E = sqrt(s);
  return ( 2.*ks/(M_PI*E)*log((E+2.*ks)/(2.*m_m_pi)) );
}

double GS_Resonance::f(const double & s) { 
  return ( m_OSwidth*m_OSmass2/pow(m_lambda,3) *
	   ( sqr(Lambda(s,m_decmasses2[0],m_decmasses2[1])) * (h(s)-m_h0) +
	     (m_OSmass2-s) * sqr(m_lambda) * m_dh0 )  );
}


 
