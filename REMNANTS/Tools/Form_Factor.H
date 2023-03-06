#ifndef SHERPA_FORM_FACTOR_H
#define SHERPA_FORM_FACTOR_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"

namespace REMNANTS {
  struct matter_form {
    enum class code {
      Single_Gaussian = 1,
      Double_Gaussian = 2,
      unknown         = -1
    };
  };
  std::ostream &operator<<(std::ostream&, const matter_form::code&);
  std::istream &operator>>(std::istream&, matter_form::code&);

  
  class Form_Factor {
  private:
    matter_form::code m_form;
    double m_fraction1, m_radius1, m_radius2;
  public:
    Form_Factor(const ATOOLS::Flavour & flav);
    ~Form_Factor() {};

    inline void
    FillParameters(matter_form::code & form, double & fraction1, double & radius1, double & radius2) {
      form = m_form; fraction1 = m_fraction1; m_radius1 = radius1; m_radius2 = radius2;
    }
    ATOOLS::Vec4D operator()();
  };
}

#endif
