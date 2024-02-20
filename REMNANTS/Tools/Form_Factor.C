#include "REMNANTS/Tools/Form_Factor.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;

Form_Factor::Form_Factor(const Flavour & flav) :
  m_flav(flav), m_form(mform::code::Single_Gaussian),
  m_fraction1(1.), m_radius1(1.), m_radius2(0.)
{
  m_form      = rempars->GetMatterForm(m_flav);
  m_radius1   = (*rempars)(flav,"MATTER_RADIUS1"); 
  m_radius2   = (m_form==mform::Single_Gaussian) ? 0. : (*rempars)(flav,"MATTER_RADIUS2"); 
  m_fraction1 = (m_form==mform::Single_Gaussian) ? 1. : (*rempars)(flav,"MATTER_FRACTION1");
  msg_Out()<<METHOD<<"("<<m_form<<": "
	   <<m_radius1<<", "<<m_radius2<<", "<<m_fraction1<<").\n";
}

Vec4D Form_Factor::operator()() {
  // Generate a position distributed according to the form-factor
  double radius = (m_form==mform::code::Double_Gaussian &&
		   ran->Get()<=m_fraction1) ? m_radius1 : m_radius2;
  double x1 = ran->GetGaussian(), x2 = ran->GetGaussian();
  return Vec4D(0.,radius*x1,radius*x2,0.);
}
