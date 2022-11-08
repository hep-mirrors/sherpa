#ifndef SHERPA_FORM_FACTOR_H
#define SHERPA_FORM_FACTOR_H

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"

namespace REMNANTS {
struct overlap_form {
  enum class code {
    Single_Gaussian = 1,
    Double_Gaussian = 2,
    unknown         = -1
  };
};


class Form_Factor {
private:
  overlap_form::code m_overlapform;
  double m_radius1, m_radius2, m_radius3, m_fraction1;
public:
  Form_Factor();
  ~Form_Factor() {};

  ATOOLS::Vec4D operator()();
};
}

#endif
