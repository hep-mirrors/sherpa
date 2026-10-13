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
  if (std::abs(s)<1.e-12) return Complex(1.,0.);
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
// Compound propagators
//
///////////////////////////////////////////////////////////////////////////

Conformal_Polynomial::Conformal_Polynomial(const std::vector<double> & c,
					   const double & s_inel,
					   const double & s1) :
  Propagator_Base(NULL,resonance_type::bespoke),
  m_c(c), m_sinel(s_inel), m_s1(s1)
{
  if (m_sinel<=m_s1)
    THROW(fatal_error,"Conformal_Polynomial needs s_inel > s_1 so that the "
	  "outer root sqrt(s_inel - s_1) is real.");
  if (m_sinel<=0.)
    THROW(fatal_error,"Conformal_Polynomial needs s_inel > 0 so that the "
	  "subtraction point z(0) is real.");
  m_outer = sqrt(m_sinel-m_s1);
  const double in0 = sqrt(m_sinel);
  m_z0 = (m_outer-in0)/(m_outer+in0);
  // Not a resonance: there is no pole mass to report, and leaving m_M at a
  // stale value would make it look like one to anything that inspects it.
  m_M = 0.; m_M2 = 0.;
}

Complex Conformal_Polynomial::Z(const double & s) const {
  // sqrt(s_inel - s), with the physical-sheet branch above s_inel:
  // s -> s + i eps sends sqrt(s_inel - s - i eps) to -i sqrt(s - s_inel),
  // which puts Im z > 0 just above the inelastic threshold.  The other branch
  // is the complex conjugate and flips the sign of Im C_p.
  const Complex inner = (s<=m_sinel) ? Complex(sqrt(m_sinel-s),0.)
                                     : Complex(0.,-sqrt(s-m_sinel));
  return (m_outer-inner)/(m_outer+inner);
}

const Complex Conformal_Polynomial::operator()(const double & s) {
  const Complex z = Z(s);
  Complex zi(1.,0.), res(0.,0.);
  double  z0i(1.);
  // The sum runs from i = 1, so step both powers before using them.
  for (size_t i(0);i<m_c.size();++i) {
    zi *= z; z0i *= m_z0;
    res += m_c[i]*(zi-z0i);
  }
  return res;
}

const Complex Conformal_Polynomial::Normalised(const double & s) {
  return (*this)(s);
}

const double Conformal_Polynomial::Normalised2(const double & s) {
  return norm((*this)(s));
}

double Conformal_Polynomial::AtInfinity() const {
  double res(0.), z0i(1.), zi(1.);
  for (size_t i(0);i<m_c.size();++i) {
    zi *= -1.; z0i *= m_z0;
    res += m_c[i]*(zi-z0i);
  }
  return res;
}

void Conformal_Polynomial::SetCoefficient(const size_t & i,const double & c) {
  if (i>=m_c.size())
    THROW(fatal_error,"Conformal_Polynomial::SetCoefficient: index "+
	  std::to_string(i)+" but only "+std::to_string(m_c.size())+
	  " coefficients.");
  m_c[i] = c;
}

Subtracted_Propagator::Subtracted_Propagator(const Complex & c0) :
  Propagator_Base(NULL), m_c0(c0) {}

Subtracted_Propagator::~Subtracted_Propagator() {
  while (!m_terms.empty()) { delete m_terms.back().p_prop; m_terms.pop_back(); }
}

void Subtracted_Propagator::Add(Propagator_Base * prop,const Complex & weight) {
  if (prop==NULL) return;
  // Evaluate the subtraction point now, once.  It is 1 for a Breit-Wigner but
  // not for Gounaris-Sakurai, so it must be measured rather than assumed.
  Term t; t.p_prop = prop; t.m_weight = weight; t.m_atzero = (*prop)(0.);
  m_terms.push_back(t);
}

const Complex Subtracted_Propagator::operator()(const double & s) {
  Complex res(m_c0);
  for (size_t i(0);i<m_terms.size();++i)
    res += m_terms[i].m_weight*((*m_terms[i].p_prop)(s)-m_terms[i].m_atzero);
  return res;
}

const Complex Subtracted_Propagator::Normalised(const double & s) {
  return (*this)(s);
}

const double Subtracted_Propagator::Normalised2(const double & s) {
  return norm((*this)(s));
}

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
  if (m_normalise && std::abs(m_norm)<1.e-14) return result;
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
