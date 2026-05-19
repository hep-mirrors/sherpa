#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;


Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(matter_form::single_gaussian),
  m_radius1(1.), m_radius2(0.), m_fraction1(1.), m_softexp(0.)
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
  msg_Out()<<METHOD<<"("<<m_flav<<"): "
	   <<"R = "<<m_radius1<<" mm, "
	   <<"expo = "<<m_softexp<<".\n";
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
  /////////////////////////////////////////////////////////////////////////////
  // Calculate effective transverse radius for the form factor, given in femtometers
  /////////////////////////////////////////////////////////////////////////////
  double radius = m_radius1;

  /////////////////////////////////////////////////////////////////////////////
  // Model 1: Double Gaussian
  /////////////////////////////////////////////////////////////////////////////
  if (m_form == matter_form::double_gaussian && ran->Get() >= m_fraction1) {
    radius = m_radius2;
  }
  /////////////////////////////////////////////////////////////////////////////
  // Model 2: x-dependent Gaussian
  // R(x) = R_0 x^(-α) reflects increasing transverse size at small x
  /////////////////////////////////////////////////////////////////////////////
  else if (m_form == matter_form::x_dependent_gaussian) {
    if (x < 0. || x > 1.) {
      radius = m_radius1;
    } else {
      radius = m_radius1 * pow(1./std::max(x, 1e-6), m_softexp);
    }
  }

  return radius;
}
