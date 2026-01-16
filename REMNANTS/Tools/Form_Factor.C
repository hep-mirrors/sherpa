#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm>

using namespace REMNANTS;
using namespace ATOOLS;


Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(matter_form::single_gaussian),
  m_GeV_fm(rpa->hBar()*rpa->c()*1.e12), 
  m_radius1(1.), m_radius2(0.), m_fraction1(1.), m_softexp(0.),
  m_variation_size(1)
{
  Initialise();
}

void Form_Factor::Initialise()
{
  /////////////////////////////////////////////////////////////////////////////
  // Radii given in fm - must be translated into 1/GeV
  /////////////////////////////////////////////////////////////////////////////  
  m_form        = rempars->Matter_Form(m_flav);
  m_radius1     = rempars->Get(m_flav,"MATTER_RADIUS_1")/m_GeV_fm;
  m_fraction1   = 1.;
  if (m_form==matter_form::double_gaussian) {
    m_radius2   = rempars->Get(m_flav,"MATTER_RADIUS_2")/m_GeV_fm;
    m_fraction1 = rempars->Get(m_flav,"MATTER_FRACTION_1");
  }
  else if (m_form==matter_form::x_dependent_gaussian) {
    m_softexp  = rempars->Get(m_flav,"SOFT_EXPONENT");
  }

  m_fraction_variations = rempars->GetVariationVector(m_flav,"MATTER_FRACTION_1");
  m_radius1_variations  = rempars->GetVariationVector(m_flav,"MATTER_RADIUS_1");
  m_radius2_variations  = rempars->GetVariationVector(m_flav,"MATTER_RADIUS_2");

  m_variation_size = std::max({m_fraction_variations.size(),
                               m_radius1_variations.size(),
                               m_radius2_variations.size(),
                               size_t(1)});
  const double frac_nom   = m_fraction_variations.empty() ? m_fraction1 : m_fraction_variations.front();
  const double radius1_nom = m_radius1_variations.empty() ? m_radius1*m_GeV_fm : m_radius1_variations.front();
  const double radius2_nom = m_radius2_variations.empty() ? m_radius2*m_GeV_fm : m_radius2_variations.front();

  m_fraction_variations.resize(m_variation_size, frac_nom);
  m_radius1_variations.resize(m_variation_size, radius1_nom);
  m_radius2_variations.resize(m_variation_size, radius2_nom);

  msg_Out()<<METHOD<<"("<<m_flav<<"): "
	   <<"R = "<<m_radius1<<" 1/GeV = "<<(m_radius1*m_GeV_fm)<<" fm, "
	   <<"expo = "<<m_softexp<<".\n";
}

double Form_Factor::Fraction1At(size_t i) const {
  if (m_form!=matter_form::double_gaussian) return 1.0;
  if (i>=m_fraction_variations.size()) return m_fraction1;
  return m_fraction_variations[i];
}

double Form_Factor::Radius1At(size_t i) const {
  if (i>=m_radius1_variations.size()) return m_radius1;
  return m_radius1_variations[i]/m_GeV_fm;
}

double Form_Factor::Radius2At(size_t i) const {
  if (i>=m_radius2_variations.size()) return m_radius2;
  return m_radius2_variations[i]/m_GeV_fm;
}

double Form_Factor::B(const double & x, const double & Q2) {
  /////////////////////////////////////////////////////////////////////////////
  // Genuinely we have different forms or combinations of Gaussians, therefore
  // an integral d^2b F(b) -> db^2 F(b^2) -> db^2 exp(-b^2/R^2)
  /////////////////////////////////////////////////////////////////////////////
  return Radius(x, Q2) * sqrt(-log(ran->Get()));
}

Vec4D Form_Factor::operator()(const double & x, const double & Q2) {
  /////////////////////////////////////////////////////////////////////////////
  // Generate a position distributed according to the form-factor.
  /////////////////////////////////////////////////////////////////////////////
  double phi = 2.*M_PI*ran->Get();
  return B(x,Q2)*Vec4D(0.,cos(phi),sin(phi),0.);
}

const double Form_Factor::Radius(const double & x,const double & Q2) const {
  /////////////////////////////////////////////////////////////////////////////
  // Default: assume single Gaussian, then the radius is just the radius
  /////////////////////////////////////////////////////////////////////////////
  double radius = m_radius1;
  /////////////////////////////////////////////////////////////////////////////
  // If double Gaussian, account for the "other" Gaussian, with probability
  // of 1 - (matter fraction)
  /////////////////////////////////////////////////////////////////////////////
  if (m_form==matter_form::double_gaussian) {
    if (m_form==matter_form::double_gaussian && ran->Get()>=m_fraction1)
    radius = m_radius2;
  }
  /////////////////////////////////////////////////////////////////////////////
  // If dynamic form factor, take into account the parton cloud
  /////////////////////////////////////////////////////////////////////////////
  else if (m_form==matter_form::x_dependent_gaussian) {
    if (x<0. || x>1.) radius = m_radius1;
    else radius = m_radius1 * pow(1./x,m_softexp);
  }
  return radius;
}
