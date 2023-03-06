#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Exception.H"

std::ostream& REMNANTS::operator<<(std::ostream& s, const matter_form::code& f)
{
  switch (f) {
    case matter_form::code::Single_Gaussian: return s << "Single_Gaussian";
    case matter_form::code::Double_Gaussian: return s << "Double_Gaussian";
  }
  return s;
}

std::istream& REMNANTS::operator>>(std::istream& s, matter_form::code& f)
{
  std::string tag;
  s >> tag;
  if (tag == "Single_Gaussian")
    f = matter_form::code::Single_Gaussian;
  else if (tag == "Double_Gaussian")
    f = matter_form::code::Double_Gaussian;
  else
    THROW(fatal_error, "Unknown matter form \"" + tag + "\"");
  return s;
}


/////////////////////////////////////////////////////////////////////////////////
// TODO: construct defaults from a flavour table and create a meaningful input
// strategy of the parameters below.
/////////////////////////////////////////////////////////////////////////////////

REMNANTS::Form_Factor::Form_Factor(const ATOOLS::Flavour & flav) :
  m_form(matter_form::code::Single_Gaussian),
  m_fraction1(1.), m_radius1(1.), m_radius2(0.)
{
  if (flav.IsNucleon()) {
    /////////////////////////////////////////////////////////////////////////////
    // we use the mean charge radius of the proton here - may need to revise
    /////////////////////////////////////////////////////////////////////////////
    m_form      = matter_form::code::Single_Gaussian;
    m_fraction1 = 1.;
    m_radius1   = 0.86;
    m_radius2   = 0.00;
  }
  else if (flav.IsPhoton()) {
    /////////////////////////////////////////////////////////////////////////////
    // we use the mean charge radius of the rho here - may need to revise
    /////////////////////////////////////////////////////////////////////////////
    m_form      = matter_form::code::Single_Gaussian;
    m_fraction1 = 1.;
    m_radius1   = 0.75;
    m_radius2   = 0.00;
  }
}

ATOOLS::Vec4D REMNANTS::Form_Factor::operator()() {
  ///////////////////////////////////////////////////////////////////////////////
  // Generate a position distributed according to the form-factor 
  ///////////////////////////////////////////////////////////////////////////////
  double radius = (m_form==matter_form::code::Double_Gaussian &&
		   ATOOLS::ran->Get()>m_fraction1) ? m_radius2 : m_radius1;
  double x1 = ATOOLS::ran->GetGaussian(), x2 = ATOOLS::ran->GetGaussian();
  return ATOOLS::Vec4D(0.,radius*x1,radius*x2,0.);
}
