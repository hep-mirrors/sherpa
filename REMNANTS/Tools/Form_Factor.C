#include "REMNANTS/Tools/Form_Factor.H"

// todo here: make sure we can also have double Gaussian and parameters
// from scoped settings.

REMNANTS::Form_Factor::Form_Factor() :
                             m_overlapform(overlap_form::Single_Gaussian),
                             m_fraction1(1.), m_radius1(1.)
{}

ATOOLS::Vec4D REMNANTS::Form_Factor::operator()() {
  double radius = m_radius1;
  if (m_overlapform==overlap_form::Double_Gaussian) {
    double rand = ATOOLS::ran->Get()-ATOOLS::sqr(m_fraction1);
    if (rand>=0.) {
      if ((rand-=ATOOLS::sqr(1-m_fraction1))<=0.) radius = m_radius2;
      else                                radius = m_radius3;
    }
  }
  double x1 = ATOOLS::ran->GetGaussian(), x2 = ATOOLS::ran->GetGaussian();
  return ATOOLS::Vec4D(0.,radius*x1,radius*x2,0.);
}