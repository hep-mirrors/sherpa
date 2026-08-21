#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"
#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Message.H"
#include <cmath>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Propagator_Base::Propagator_Base(Total_Width_Base * width,
				 const resonance_type & type) :
  p_width(width), m_type(type), m_M(0.), m_M2(m_M*m_M) {
  if (p_width!=NULL) { m_M = p_width->Flav().Mass(true); m_M2 = m_M*m_M; } 
}

double Propagator_Base::OnShellWidth() const {
  return (p_width!=NULL ? (*p_width)(m_M2) : 0.);
}

ATOOLS::Flavour Propagator_Base::Flav() const {
  return (p_width!=NULL ? p_width->Flav() : ATOOLS::Flavour(kf_none));
}


///////////////////////////////////////////////////////////////////////////
//
// Simple Breit Wigner (m_type==fixed/running), and Gounaris-Sakurai
// (m_type==GS).
//
// The GS form replaces the constant numerator M^2 and the plain
// s-independent M^2-s denominator term by:
//   BW^GS(s) = (M^2 + d*M*Gamma(s)) /
//              (M^2 - s + f(s) - i*sqrt(s)*Gamma(s))
// with f(s), d and the k(s), h(s) helper functions as defined in
// Gounaris & Sakurai, PRL 21 (1968) 244, using the same conventions as
// arXiv:1509.09140, Eq.(2.2) (BW^GS). Gamma(s) is taken from the same
// (running) Total_Width_Base used for the plain Breit-Wigner, i.e. we
// re-use whatever partial-width machinery already exists instead of
// re-deriving Gamma_pipi(s) here.
//
///////////////////////////////////////////////////////////////////////////

BreitWigner::BreitWigner(Total_Width_Base * width,
			 const resonance_type & type,
			 const double & mdau) :
  Propagator_Base(width,type),
  m_mdau2(mdau>0. ? mdau*mdau : sqr(ATOOLS::Flavour(kf_pi_plus).Mass(true))),
  m_GammaPole(0.), m_kPole(0.), m_hPole(0.), m_hPrimePole(0.), m_d(0.)
{
  if (m_type==resonance_type::GS) InitGS();
}

double BreitWigner::k(const double & s) const {
  double arg = s/4.-m_mdau2;
  return (arg>0. ? sqrt(arg) : 0.);
}

double BreitWigner::h(const double & s) const {
  double ks = k(s);
  if (s<=4.*m_mdau2 || ks<=0.) return 0.;
  return (2./M_PI)*(ks/sqrt(s))*log((sqrt(s)+2.*ks)/(2.*sqrt(m_mdau2)));
}

void BreitWigner::InitGS() {
  if (p_width==NULL) return;
  m_GammaPole = (*p_width)(m_M2);
  m_kPole     = k(m_M2);
  if (m_kPole<=0.) {
    msg_Error()<<"Error in "<<METHOD<<": pole mass below 2*mdau threshold, "
	       <<"Gounaris-Sakurai parametrization undefined. "
	       <<"Falling back to f(s)=0, d=0 (i.e. plain running-width BW).\n";
    return;
  }
  m_hPole      = h(m_M2);
  m_hPrimePole = m_hPole*(1./(8.*sqr(m_kPole))-1./(2.*m_M2)) +
                 1./(2.*M_PI*m_M2);
  m_d          = (3./M_PI)*(m_mdau2/sqr(m_kPole)) *
                 log((m_M+2.*m_kPole)/(2.*sqrt(m_mdau2))) +
                 m_M/(2.*M_PI*m_kPole) -
                 m_mdau2*m_M/(M_PI*pow(m_kPole,3));
}

const Complex BreitWigner::ValueGS(const double & s) {
  double Gs = (*p_width)(s);
  double fs = 0.;
  if (m_kPole>0.) {
    double ks = k(s);
    fs = m_GammaPole*m_M2/pow(m_kPole,3) *
         ( sqr(ks)*(h(s)-m_hPole) + (m_M2-s)*sqr(m_kPole)*m_hPrimePole );
  }
  Complex num(m_M2 + m_d*m_M*Gs, 0.);
  Complex den(m_M2 - s + fs, -sqrt(s)*Gs);
  return num/den;
}

const Complex BreitWigner::operator()(const double & s) {
  if (m_type==resonance_type::GS) return ValueGS(s);
  return m_M2/Complex(m_M2-s,-sqrt(s)*(*p_width)(s));
}

const double BreitWigner::Normalised2(const double & s) {
  if (m_type==resonance_type::GS) return norm(ValueGS(s));
  double Gamma = (*p_width)(s), MG2 = s*sqr(Gamma);
  return m_M2/(sqr(s-m_M2)+MG2);
}

const Complex BreitWigner::Normalised(const double & s) {
  if (m_type==resonance_type::GS) return ValueGS(s);
  return m_M2/Complex(m_M2-s,-sqrt(s)*(*p_width)(s));
}


///////////////////////////////////////////////////////////////////////////
//
// RChL_BW: see the class comment in Propagator.H for the convention.
//
///////////////////////////////////////////////////////////////////////////

const Complex RChL_BW::operator()(const double & s) {
  return 1./Complex(s-m_M2,-m_M*(*p_width)(s));
}

const Complex RChL_BW::Normalised(const double & s) {
  return (*this)(s);
}

const double RChL_BW::Normalised2(const double & s) {
  return norm((*this)(s));
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

///////////////////////////////////////////////////////////////////////////
//
// Diagnostic dump (request #1). Handles three cases: a Summed_Propagator
// (prints each constituent's mass/on-shell width/weight), a single
// Propagator_Base (prints just that one), or NULL/unrecognised (prints
// that nothing was constructed) - all with the same one-line-per-entry
// format so grep'ing the log for a given channel is straightforward.
//
///////////////////////////////////////////////////////////////////////////

void METOOLS::DumpPropagatorStructure(const std::string & label,
				       const int & ffmodel_id,
				       Propagator_Base * props) {
  msg_Out()<<"### Propagator structure for \""<<label<<"\" "
	   <<"(FORM_FACTOR = "<<ffmodel_id<<"):\n";
  if (props==NULL) {
    msg_Out()<<"###   <none constructed - falls back to a constant "
	     <<"form factor>\n";
    return;
  }
  Summed_Propagator * summed = dynamic_cast<Summed_Propagator *>(props);
  if (summed!=NULL) {
    for (map<Propagator_Base *,Complex>::iterator pit=summed->GetAll().begin();
	 pit!=summed->GetAll().end();pit++) {
      Propagator_Base * p = pit->first;
      msg_Out()<<"###   "<<p->Flav()<<":  M = "<<p->Mass()<<" GeV,  "
	       <<"Gamma(M^2) = "<<p->OnShellWidth()<<" GeV,  "
	       <<"weight = "<<pit->second<<"\n";
    }
    return;
  }
  msg_Out()<<"###   "<<props->Flav()<<":  M = "<<props->Mass()<<" GeV,  "
	   <<"Gamma(M^2) = "<<props->OnShellWidth()<<" GeV\n";
}
