#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;


Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(matter_form::single_gaussian),
  m_radius1(1.), m_radius2(0.), m_fraction1(1.), m_softexp(0.),
  m_n_matter_form_variations(1)
{
  Initialise();
}

void Form_Factor::Initialise()
{
  /////////////////////////////////////////////////////////////////////////////
  // Radii given in fm - must be translated into mm
  /////////////////////////////////////////////////////////////////////////////
  m_form        = rempars->Matter_Form(m_flav);
  m_radius1     = rempars->Get(m_flav,"MATTER_RADIUS_1");
  m_fraction1   = 1.;
  if (m_form==matter_form::double_gaussian) {
    m_radius2   = rempars->Get(m_flav,"MATTER_RADIUS_2");
    m_fraction1 = rempars->Get(m_flav,"MATTER_FRACTION_1");
  }
  else if (m_form==matter_form::x_dependent_gaussian) {
    m_softexp  = rempars->Get(m_flav,"SOFT_EXPONENT");
  }

  m_fraction_variations = rempars->GetVariationVector(m_flav,"MATTER_FRACTION_1");
  m_radius1_variations  = rempars->GetVariationVector(m_flav,"MATTER_RADIUS_1");
  m_radius2_variations  = rempars->GetVariationVector(m_flav,"MATTER_RADIUS_2");
  m_softexp_variations  = rempars->GetVariationVector(m_flav,"SOFT_EXPONENT");

  m_n_matter_form_variations = std::max({m_fraction_variations.size(),
                               m_radius1_variations.size(),
                               m_radius2_variations.size(),
                               m_softexp_variations.size(),
                               size_t(1)});
  if (m_n_matter_form_variations > 1) {
    const double frac_nom    = m_fraction_variations.empty() ? m_fraction1 : m_fraction_variations.front();
    const double radius1_nom = m_radius1_variations.empty() ? m_radius1 : m_radius1_variations.front();
    const double radius2_nom = m_radius2_variations.empty() ? m_radius2 : m_radius2_variations.front();
    const double softexp_nom = m_softexp_variations.empty() ? m_softexp : m_softexp_variations.front();

    m_fraction_variations.resize(m_n_matter_form_variations, frac_nom);
    m_radius1_variations.resize(m_n_matter_form_variations, radius1_nom);
    m_radius2_variations.resize(m_n_matter_form_variations, radius2_nom);
    m_softexp_variations.resize(m_n_matter_form_variations, softexp_nom);
  }

  msg_Out()<<METHOD<<"("<<m_flav<<"): "
	   <<"R = "<<m_radius1<<" mm, "
	   <<"expo = "<<m_softexp<<".\n";
}

double Form_Factor::Fraction1At(size_t ivar) const {
  if (m_form!=matter_form::double_gaussian) return 1.0;
  if (m_n_matter_form_variations<=1 || ivar>=m_n_matter_form_variations) return m_fraction1;
  return m_fraction_variations[ivar];
}

double Form_Factor::Radius1At(size_t ivar) const {
  if (m_n_matter_form_variations<=1 || ivar>=m_n_matter_form_variations) return m_radius1;
  return m_radius1_variations[ivar];
}

double Form_Factor::Radius2At(size_t ivar) const {
  if (m_form!=matter_form::double_gaussian) return 0.0;
  if (m_n_matter_form_variations<=1 || ivar>=m_n_matter_form_variations) return m_radius2;
  return m_radius2_variations[ivar];
}

double Form_Factor::SoftExponentAt(size_t ivar) const {
  if (m_form!=matter_form::x_dependent_gaussian) return 0.0;
  if (m_n_matter_form_variations<=1 || ivar>=m_n_matter_form_variations) return m_softexp;
  return m_softexp_variations[ivar];
}

double Form_Factor::B(const double & x, const double & Q2) {
  /////////////////////////////////////////////////////////////////////////////
  // Sample radial distance from 2D Gaussian form factor F(b) = exp(-b^2/R^2)
  // In polar coordinates: P(b)db = b exp(-b^2/R^2)db
  // Inverse transform sampling gives: b = R sqrt(-ln(u)) where u ~ Uniform(0,1)
  // Returns impact parameter b in units of millimeter.
  /////////////////////////////////////////////////////////////////////////////
  return Radius(x, Q2) * sqrt(-log(ran->Get())) * 1.e-12;
}

Vec4D Form_Factor::operator()(const double & x, const double & Q2) {
  /////////////////////////////////////////////////////////////////////////////
  // Generate a 2D position distributed according to the form factor.
  // Position is in the transverse (x-y) plane relative to the parton center.
  // Azimuthal angle \phi is uniformly distributed in [0, 2\pi).
  // Returns position in units of millimeter.
  /////////////////////////////////////////////////////////////////////////////
  double phi = 2. * M_PI * ran->Get();
  return B(x, Q2) * Vec4D(0., std::cos(phi), std::sin(phi), 0.);
}

const double Form_Factor::Radius(const double & x, const double & Q2) const {
  return Radius(x, Q2, 0);
}

const double Form_Factor::Radius(const double & x, const double & Q2,
                                 size_t ivar) const {
  /////////////////////////////////////////////////////////////////////////////
  // Calculate effective transverse radius for the form factor, given in femtometers
  /////////////////////////////////////////////////////////////////////////////
  const double radius1 = Radius1At(ivar);
  double radius = radius1;

  /////////////////////////////////////////////////////////////////////////////
  // Model 1: Double Gaussian
  /////////////////////////////////////////////////////////////////////////////
  if (m_form == matter_form::double_gaussian &&
      ran->Get() >= Fraction1At(ivar)) {
    radius = Radius2At(ivar);
  }
  /////////////////////////////////////////////////////////////////////////////
  // Model 2: x-dependent Gaussian
  // R(x) = R_0 x^(-α) reflects increasing transverse size at small x
  /////////////////////////////////////////////////////////////////////////////
  else if (m_form == matter_form::x_dependent_gaussian) {
    if (x < 0. || x > 1.) {
      radius = radius1;
    } else {
      radius = radius1 * pow(1./std::max(x, 1e-6), SoftExponentAt(ivar));
    }
  }

  return radius;
}
