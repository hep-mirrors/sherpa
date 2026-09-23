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
// Complex-pole propagator. See Propagator.H for why a Breit-Wigner is
// not an acceptable substitute for the very broad scalars.
//
// PDG T-matrix poles, sqrt(s_pole) = M - i Gamma/2:
//   f_0(500)/sigma     M = 449 (+22-16) MeV,  Gamma/2 = 275 +- 12 MeV
//   K*_0(700)/kappa    M = 680 +- 50 MeV,     Gamma/2 = 300 +- 40 MeV
// i.e. Gamma > M for the sigma. Note that Gamma here is the FULL pole
// width, so the imaginary part of sqrt(s_pole) is Gamma/2.
//
///////////////////////////////////////////////////////////////////////////

Complex_Pole::Complex_Pole(const double & Mpole,const double & Gpole,
                           const ATOOLS::Flavour & flav) :
  Propagator_Base(NULL,resonance_type::complex_pole),
  m_Mpole(Mpole), m_Gpole(Gpole), m_flav(flav) {
  const Complex sqrts(m_Mpole,-0.5*m_Gpole);
  m_spole = sqrts*sqrts;
  // m_M/m_M2 are left as the REAL part of the pole position purely so
  // that Mass() reports something sensible in a dump; nothing in the
  // propagator itself uses them.
  m_M  = m_Mpole;
  m_M2 = m_M*m_M;
}

const Complex Complex_Pole::operator()(const double & s) {
  return m_spole/(m_spole-s);
}

const Complex Complex_Pole::Normalised(const double & s) {
  return (*this)(s);
}

const double Complex_Pole::Normalised2(const double & s) {
  const Complex v = (*this)(s);
  return (v*conj(v)).real();
}


///////////////////////////////////////////////////////////////////////////
//
// Flatte. See Propagator.H for why a Breit-Wigner will not do for a
// resonance sitting on a threshold.
//
// Standard parameter sets, couplings in GeV^2:
//   f_0(980), BES:            M = 0.965, g_pipi = 0.165, g_KK = 0.695
//                             (g_KK/g_pipi = 4.21)
//   a_0(980), Crystal Barrel: M = 0.999, g_etapi = 0.221, g_KK = 0.256
// Note that these M values are NOT the PDG "mass": for a Flatte the
// parameters are correlated and only meaningful as a set.
//
///////////////////////////////////////////////////////////////////////////

Flatte::Flatte(const double & M,const double & g1,const double & g2,
               const double & ma1,const double & mb1,
               const double & ma2,const double & mb2,
               const ATOOLS::Flavour & flav) :
  Propagator_Base(NULL,resonance_type::flatte),
  m_g1(g1), m_g2(g2),
  m_ma1(ma1), m_mb1(mb1), m_ma2(ma2), m_mb2(mb2), m_flav(flav) {
  m_M = M; m_M2 = M*M;
}

Complex Flatte::Rho(const double & s,const double & ma,
                    const double & mb) const {
  if (s<=0.) return Complex(0.,0.);
  const double a = 1.-sqr(ma+mb)/s, b = 1.-sqr(ma-mb)/s;
  const double r2 = a*b;
  // Above threshold rho is real; below it the square root is imaginary
  // and the channel feeds the REAL part of the denominator instead.
  if (r2>=0.) return Complex(sqrt(r2),0.);
  return Complex(0.,sqrt(-r2));
}

const Complex Flatte::operator()(const double & s) {
  const Complex den = Complex(m_M2-s,0.)
    - Complex(0.,1.)*(m_g1*Rho(s,m_ma1,m_mb1)+m_g2*Rho(s,m_ma2,m_mb2));
  if (std::abs(den)<1.e-30) return Complex(0.,0.);
  return m_M2/den;
}

const Complex Flatte::Normalised(const double & s) { return (*this)(s); }

const double Flatte::Normalised2(const double & s) {
  const Complex v = (*this)(s);
  return (v*conj(v)).real();
}

double Flatte::OnShellWidth() const {
  // Gamma = (g_1 rho_1 + g_2 rho_2)/M at s=M^2, counting only the
  // channels that are actually open there.
  const Complex r1 = Rho(m_M2,m_ma1,m_mb1), r2 = Rho(m_M2,m_ma2,m_mb2);
  return (m_g1*r1.real()+m_g2*r2.real())/m_M;
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

///////////////////////////////////////////////////////////////////////////
//
// Diagnostic dump (request #1). Recurses into nested Summed_Propagator/
// Multiplied_Propagator constituents (e.g. the "rho(770)+alpha*(rho x
// omega)" combination used for pi- pi+ pi-'s rho-omega mixing term) -
// a constituent that is ITSELF a composite propagator has no single
// Flavour/mass of its own (Propagator_Base::Flav() falls back to
// kf_none="no_particle" for it), so printing it as if it were a leaf
// resonance is misleading (looks like a missing/unregistered particle
// when the underlying physics is actually fine - confirmed by an
// explicit review: the composite's own leaf constituents, once
// recursed into, are properly registered). Handles three cases per
// level: a Summed_Propagator, a Multiplied_Propagator, or a leaf
// Propagator_Base with a real Flavour.
//
///////////////////////////////////////////////////////////////////////////

static void DumpPropagatorEntry(Propagator_Base * p, const Complex & weight,
				 const std::string & indent, bool haveWeight) {
  Summed_Propagator     * sub_s = dynamic_cast<Summed_Propagator     *>(p);
  Multiplied_Propagator  * sub_m = dynamic_cast<Multiplied_Propagator *>(p);
  if (sub_s!=NULL || sub_m!=NULL) {
    msg_Out()<<"###   "<<indent
	     <<(sub_s!=NULL ? "[nested sum]" : "[nested product]");
    if (haveWeight) msg_Out()<<", weight = "<<weight;
    msg_Out()<<"\n";
    map<Propagator_Base *,Complex> & sub =
      (sub_s!=NULL ? sub_s->GetAll() : sub_m->GetAll());
    for (map<Propagator_Base *,Complex>::iterator sit=sub.begin();
	 sit!=sub.end();sit++) {
      DumpPropagatorEntry(sit->first, sit->second, indent+"  ", true);
    }
    return;
  }
  msg_Out()<<"###   "<<indent<<p->Flav()<<":  M = "<<p->Mass()<<" GeV,  "
	   <<"Gamma(M^2) = "<<p->OnShellWidth()<<" GeV";
  if (haveWeight) msg_Out()<<",  weight = "<<weight;
  msg_Out()<<"\n";
}

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
  Summed_Propagator      * summed = dynamic_cast<Summed_Propagator     *>(props);
  Multiplied_Propagator   * mult  = dynamic_cast<Multiplied_Propagator *>(props);
  if (summed!=NULL || mult!=NULL) {
    map<Propagator_Base *,Complex> & top =
      (summed!=NULL ? summed->GetAll() : mult->GetAll());
    for (map<Propagator_Base *,Complex>::iterator pit=top.begin();
	 pit!=top.end();pit++) {
      DumpPropagatorEntry(pit->first, pit->second, "", true);
    }
    return;
  }
  DumpPropagatorEntry(props, Complex(1.,0.), "", false);
}
