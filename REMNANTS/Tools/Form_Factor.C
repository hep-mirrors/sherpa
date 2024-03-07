#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;


Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(matter_form::single_gaussian),
  m_radius1(1.), m_radius2(0.), m_fraction1(1.)
  
{
  Initialise();
}

void Form_Factor::Initialise()
{
  m_form        = rempars->Matter_Form(m_flav);
  m_radius1     = rempars->Get(m_flav,"RADIUS_1");
  if (m_form==matter_form::single_gaussian) {
    m_fraction1 = 1.;
  }
  if (m_form==matter_form::double_gaussian) {
    m_radius2   = rempars->Get(m_flav,"RADIUS_2");
    m_fraction1 = rempars->Get(m_flav,"FRACTION_1");
  }
}

Vec4D Form_Factor::operator()() {
  // Generate a position distributed according to the form-factor
  double radius = (m_form==matter_form::double_gaussian &&
		   ran->Get()<=m_fraction1) ? m_radius1 : m_radius2;
  double x1 = ran->GetGaussian(), x2 = ran->GetGaussian();
  return Vec4D(0.,radius*x1,radius*x2,0.);
}
