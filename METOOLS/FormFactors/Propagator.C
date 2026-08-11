#include "METOOLS/FormFactors/Propagator.H"
#include "METOOLS/FormFactors/Resonance_Base.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Propagator_Base::Propagator_Base(Total_Width_Base * width,
				 const resonance_type & type,
				 const double & M,const double & Gamma) :
  p_width(width), m_type(type), m_M(0.), m_M2(0.), m_Gamma(0.) {
  if (p_width!=NULL) {
    m_M     = p_width->Flav().Mass(true);
    // Width(true): otherwise anything flagged as narrow silently comes back
    // with a vanishing width, leaving an undamped pole in the propagator.
    m_Gamma = p_width->Flav().Width(true);
  }
  // Explicitly given parameters win over the particle table.
  if (M>0.)     m_M     = M;
  if (Gamma>0.) m_Gamma = Gamma;
  m_M2 = m_M*m_M;
}

const double Propagator_Base::Width(const double & s) const {
  if (m_type==resonance_type::fixed || p_width==NULL) return m_Gamma;
  return (*p_width)(s);
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
  return s/Complex(m_M2-s,-m_M*Width(s));
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
				 const resonance_type & type,
				 const double & M,const double & Gamma) :
  Propagator_Base(width,type,M,Gamma) {
  // Line_Shapes::Get() hands back a NULL pointer for resonances it does not
  // know about, so the parameters have to be supplied explicitly then.
  if (m_M<=0. || m_Gamma<=0.)
    THROW(fatal_error,"Gounaris-Sakurai propagator without mass or width - "
	  "there is no lineshape for this resonance, pass them explicitly.");
  m_mpi   = Flavour(kf_pi_plus).Mass();
  m_mpi2  = sqr(m_mpi);
  m_ppiM  = ppi(m_M2);
  m_d     = ( 3./M_PI * m_mpi2/(m_ppiM*m_ppiM) *
	      log((m_M+2.*m_ppiM)/(2.*m_mpi)) +
	      m_M/(2.*M_PI*m_ppiM) -
	      m_M*m_mpi2/(M_PI*m_ppiM*m_ppiM*m_ppiM) );
  m_hM2   = h(m_M2);
  m_dhM2  = dh(m_M2);
}

const Complex GounarisSakurai::ppi(const double & q2) const {
  return 0.5 * sqrt(Complex(q2-4.*m_mpi2,0.));
}

const Complex GounarisSakurai::h(const double & q2) const {
  Complex p = ppi(q2), q = sqrt(Complex(q2,0.));
  return 2./M_PI * p/q * log((q+2.*p)/(2.*m_mpi));
}
const Complex GounarisSakurai::dh(const double & q2) const {
  Complex p = ppi(q2);
  return h(q2)/8. * (1./(p*p)-4./q2) + 1./(2.*M_PI*q2);
}
const Complex GounarisSakurai::f(const double & q2) const {
  Complex p = ppi(q2);
  return ( m_Gamma*m_M2/(m_ppiM*m_ppiM*m_ppiM) *
	   ( p*p*(h(q2)-m_hM2) + m_ppiM*m_ppiM*(m_M2-q2)*m_dhM2 ) );
}
const Complex GounarisSakurai::GammaGS(const double & q2) const {
  Complex r = ppi(q2)/m_ppiM;
  return m_Gamma * m_M/sqrt(q2) * r*r*r;
}

const Complex GounarisSakurai::operator()(const double & s) {
  // The d and f terms above are derived for a pure pi pi width, so the
  // imaginary part has to use the same one - mixing in a multi-channel
  // running width from a lineshape breaks the F(0) = 1 normalisation.
  Complex width = (m_type==resonance_type::GS)?
    GammaGS(s) : Complex(Width(s),0.);
  return ( (m_M2+m_d*m_M*m_Gamma)/
	   (Complex(m_M2-s,0.)+f(s)-Complex(0.,1.)*m_M*width) );
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

Summed_Propagator::Summed_Propagator(const bool & normalise) :
  Propagator_Base(NULL),
  m_norm(Complex(0.,0.)), m_normalise(normalise) {}

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
  if (m_props.size()==0) return 1.0;
  Complex result(0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    result += pit->second*(*pit->first)(s);
  }
  return m_normalise?result/m_norm:result;
}

const Complex Summed_Propagator::Normalised(const double & s) {
  Complex ampl(0.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl += pit->second*pit->first->Normalised(s);
  }
  return m_normalise?ampl/m_norm:ampl;
}

const double Summed_Propagator::Normalised2(const double & s) {
  return norm(Normalised(s));
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
  Complex ampl(1.,0.);
  for (map<Propagator_Base *,Complex>::iterator pit=m_props.begin();
       pit!=m_props.end();pit++) {
    ampl *= pit->second*pit->first->Normalised(s);
  }
  return ampl/m_norm;
}

const double Multiplied_Propagator::Normalised2(const double & s) {
  return norm(Normalised(s));
}
