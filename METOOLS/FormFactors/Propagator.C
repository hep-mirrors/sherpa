#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Resonance_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Propagator_Base::Propagator_Base(Total_Width_Base * width,
				 const resonance_type & type) :
  p_width(width), m_type(type), m_M(0.), m_M2(m_M*m_M) {
  if (p_width!=NULL) { m_M = p_width->Flav().Mass(true); m_M2 = m_M*m_M; } 
}


///////////////////////////////////////////////////////////////////////////
//
// Simple Breit Wigner
//
///////////////////////////////////////////////////////////////////////////

const Complex BreitWigner::operator()(const double & s) {
  return m_M2/Complex(m_M2-s,-sqrt(s)*(*p_width)(s));
}

const double BreitWigner::Normalised2(const double & s) {
  double Gamma = (*p_width)(s), MG2 = s*sqr(Gamma);
  return m_M2/(sqr(s-m_M2)+MG2);
}

const Complex BreitWigner::Normalised(const double & s) {
  return m_M2/Complex(m_M2-s,-m_M*(*p_width)(s));
}

///////////////////////////////////////////////////////////////////////////
//
// Weighted Breit Wigner
//
///////////////////////////////////////////////////////////////////////////

const Complex WeightedBreitWigner::operator()(const double & s) {
  return s/Complex(m_M2-s,-sqrt(s)*(*p_width)(s));
}

const Complex WeightedBreitWigner::Normalised(const double & s) {
  THROW(fatal_error,"Normalised not implemented for Weighted BW form.");
}

const double  WeightedBreitWigner::Normalised2(const double & s) {
  THROW(fatal_error,"Normalised2 not implemented for Weighted BW form.");
}

///////////////////////////////////////////////////////////////////////////
//
// Gounaris Sakurai
//
///////////////////////////////////////////////////////////////////////////

GounarisSakurai::GounarisSakurai(Total_Width_Base * width,
				 const resonance_type & type) :
  Propagator_Base(width,type) {
  msg_Out()<<METHOD<<" width = ["<<p_width<<"]\n";
  m_Gamma = p_width->Flav().Width();
  m_mpi   = Flavour(kf_pi_plus).Mass();
  m_mpi2  = sqr(m_mpi);
  m_ppi2  = (m_M2-4.*m_mpi2)/4.;
  m_ppi   = sqrt(m_ppi2);
  m_d     = ( 3./M_PI * m_mpi2/m_ppi2 * log((m_M+2.*m_ppi)/(2.*m_mpi)) +
	      m_M/(2.*M_PI*m_ppi) - m_M*m_mpi2/(M_PI*pow(m_ppi,3)) );
  m_hM2   = h(m_M2);
  m_dhM2  = dh(m_M2);
}

const double GounarisSakurai::h(const double & q2) const {
  double q = sqrt(q2);
  return 2./M_PI * m_ppi/q * log((q+2.*m_ppi)/(2.*m_mpi));
}
const double GounarisSakurai::dh(const double & q2) const {
  double ppi2 = q2/4.-m_mpi2;
  return h(q2)/8. * (1./ppi2-4./q2) + 1./(2.*M_PI*q2);
}
const double GounarisSakurai::f(const double & q2) const {
  return ( m_Gamma*m_M2/pow(m_ppi,3.) *
	   ( (q2/4.-m_mpi2)*(h(q2)-m_hM2) +
	     m_ppi2*(m_M2-q2)*dh(q2) ) );
}

const Complex GounarisSakurai::operator()(const double & s) {
  return ( (m_M2+m_d*m_M*m_Gamma)/
	   Complex(m_M2-s+f(s), -m_M*(*p_width)(s)) );
}

const Complex GounarisSakurai::Normalised(const double & s) {
  THROW(fatal_error,"Normalised not implemented for Gounaris-Sakurai form.");
}

const double  GounarisSakurai::Normalised2(const double & s) {
  THROW(fatal_error,"Normalised2 not implemented for Gounaris-Sakurai form.");
}

///////////////////////////////////////////////////////////////////////////
//
// Compound propagators, needed, e.g. for form factors 
//
///////////////////////////////////////////////////////////////////////////

Summed_Propagator::Summed_Propagator(Propagator_Base * prop) :
  Propagator_Base(NULL),
  m_norm(Complex(0.,0.)) {
  if (prop!=NULL) m_props[prop] = m_norm;
}

Summed_Propagator::~Summed_Propagator() {
  while (!m_props.empty()) {
    delete m_props.begin()->first;
    m_props.erase(m_props.begin());
  }
}
   
void Summed_Propagator::Add(Propagator_Base * prop,const Complex & weight) {
  if (prop!=NULL && m_props.find(prop)==m_props.end()) {
    m_props[prop] = weight;
    m_norm       += weight;
  }
}

const Complex Summed_Propagator::operator()(const double & s) {
  Complex result(0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    result += pit->second*(*pit->first)(s);
  }
  return result/m_norm;
}

const Complex Summed_Propagator::Normalised(const double & s) {
  Complex ampl = (0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl += pit->second*pit->first->Normalised(s);
  }
  return ampl/m_norm;
}

const double Summed_Propagator::Normalised2(const double & s) {
  Complex ampl = (0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl += pit->second*pit->first->Normalised(s);
  }
  return norm(ampl/m_norm);
}

///////////////////////////////////////////////////////////////////////////
//
// Multiplied propagators, needed, e.g. for form factors 
//
///////////////////////////////////////////////////////////////////////////

Multiplied_Propagator::Multiplied_Propagator(Propagator_Base * prop) :
  Propagator_Base(NULL),
  m_norm(Complex(1.,0.)) {
  if (prop!=NULL) m_props[prop] = m_norm;
}

Multiplied_Propagator::~Multiplied_Propagator() {
  while (!m_props.empty()) {
    delete m_props.begin()->first;
    m_props.erase(m_props.begin());
  }
}
   
void Multiplied_Propagator::Add(Propagator_Base * prop,const Complex & weight) {
  if (prop!=NULL && m_props.find(prop)==m_props.end()) {
    m_props[prop] = weight;
    m_norm       *= weight;
  }
}

const Complex Multiplied_Propagator::operator()(const double & s) {
  Complex result(1.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    result *= pit->second*(*pit->first)(s);
  }
  return result/m_norm;
}

const Complex Multiplied_Propagator::Normalised(const double & s) {
  Complex ampl = (0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl *= pit->second*pit->first->Normalised(s);
  }
  return ampl/m_norm;
}

const double Multiplied_Propagator::Normalised2(const double & s) {
  Complex ampl = (0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl *= pit->second*pit->first->Normalised(s);
  }
  return norm(ampl/m_norm);
}
