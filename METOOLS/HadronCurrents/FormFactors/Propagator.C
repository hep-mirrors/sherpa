#include "METOOLS/HadronCurrents/Tools.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"
#include "METOOLS/HadronCurrents/FormFactors/Novosibirsk4pi_GTables.H"
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

FM95_Fixed_BW::FM95_Fixed_BW(const double & mass,const double & width) :
  Propagator_Base(NULL,resonance_type::fixed), m_G(width) {
  m_M  = mass;
  m_M2 = mass*mass;
}

const Complex FM95_Fixed_BW::operator()(const double & s) {
  return Tools::BreitWignerFix(s,m_M2,m_M*m_G);
}

const Complex FM95_Fixed_BW::Normalised(const double & s) {
  return (*this)(s);
}

const double FM95_Fixed_BW::Normalised2(const double & s) {
  return norm((*this)(s));
}

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

// ====================================================================
// Novo4Pi_Propagator - see Propagator.H for why these are not BreitWigner
// with resonance_type::GS.
// ====================================================================

Novo4Pi_Propagator::Novo4Pi_Propagator(const novo4pi_resonance & res,
                                 const double & M,const double & G,
                                 const double & mdau,
                                 const Complex & z,const double & scaleA,
                                 Total_Width_Base * width,
                                 const bool & keepshape) :
  Propagator_Base(width,resonance_type::running),
  m_res(res), m_G(G), m_scaleA(scaleA), m_z(z), m_mdau(mdau),
  m_rhonorm(1.), m_g1(1.), m_keepshape(keepshape) {
  m_M = M; m_M2 = M*M;
  if (m_res==novo4pi_resonance::rho) {
    const double d0 = Novo4PiDM(0.,m_M);
    const double n  = 1.+m_G/m_M*d0;
    m_rhonorm = (std::abs(n)>1.e-12) ? 1./n : 1.;
  }
}

const Complex Novo4Pi_Propagator::operator()(const double & s) {
  // A registered line shape gives the same dimensionless form for every
  // resonance: Bondar's D with his g_R(s)/g_R(M^2) replaced by Sherpa's
  // running width.  The pole width cancels in the ratio, leaving
  // Gamma(s)/M.  For the rho this also drops the Gounaris-Sakurai
  // dispersive real part, which binp.f carries and the paper does not.
  if (p_width && m_keepshape && m_res==novo4pi_resonance::rho) {
    // Bondar's rho with only the width function replaced.
    Complex D(s-m_M2-Novo4PiDM(s,m_M)*m_G*m_M, m_M*(*p_width)(s));
    return D/m_M2*m_rhonorm;
  }
  if (p_width) return Complex(s/m_M2-1., (*p_width)(s)/m_M);
  // Only the sigma reaches here: every other resonance is constructed
  // with a registered line shape and returns above.
  switch (m_res) {
  case novo4pi_resonance::sigma: {
    // binp.f z_dsigma: equal daughter masses, so the two Kallen
    // factors coincide.
    const double d  = 1.-4.*sqr(m_mdau)/s;
    const double d0 = 1.-4.*sqr(m_mdau)/m_M2;
    const double pm  = sqrt(Max(d*d,1.e-16));
    const double pm0 = sqrt(Max(d0*d0,1.e-16));
    return Complex(s/m_M2-1., m_G/m_M*pm/pm0);
  }
  }
  return Complex(1.,0.);
}

const Complex Novo4Pi_Propagator::Normalised(const double & s) {
  // These are already dimensionless by construction; "normalised" here
  // means the same thing, so nothing further to do.
  return (*this)(s);
}

const double Novo4Pi_Propagator::Normalised2(const double & s) {
  return norm((*this)(s));
}

Two_Channel_Flatte::Two_Channel_Flatte(Total_Width_Base * width,
                                       const double & mass,
                                       const Complex & g1,const double & m11,
                                       const double & m12,
                                       const Complex & g2,const double & m21,
                                       const double & m22,
                                       const double & gamma0) :
  Propagator_Base(width,resonance_type::bespoke),
  m_g1(g1), m_m11(m11), m_m12(m12),
  m_g2(g2), m_m21(m21), m_m22(m22),
  m_gamma0(gamma0)
{
  m_M  = mass;
  m_M2 = sqr(m_M);
}

Complex Two_Channel_Flatte::PhaseSpace(const double & s,const double & m1,
                                       const double & m2) const {
  if (s==0.) return Complex(0.,0.);
  const double splus  = sqr(m1+m2);
  const double sminus = sqr(m1-m2);
  // Complex sqrt on purpose: below threshold this is imaginary, which is
  // exactly what a Flatte needs and what a running width cannot express.
  return sqrt(Complex((s-splus)*(s-sminus)/(s*s),0.));
}

const Complex Two_Channel_Flatte::operator()(const double & s) {
  const Complex sigma =
    m_g1*PhaseSpace(s,m_m11,m_m12) +
    m_g2*PhaseSpace(s,m_m21,m_m22);
  // gamma0 carries residual modes not represented by the two explicit
  // threshold channels.  The complex g_i have dimensions of GeV^2 and allow
  // phenomenological relative phases between the coupled channels.
  const Complex denom =
    Complex(m_M2-s,0.)-Complex(0.,1.)*(sqrt(s)*m_gamma0+sigma);
  return m_M2/denom;
}

const Complex Two_Channel_Flatte::Normalised(const double & s) {
  return (*this)(s);
}

const double Two_Channel_Flatte::Normalised2(const double & s) {
  return norm((*this)(s));
}
