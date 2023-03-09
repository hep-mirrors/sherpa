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
    ATOOLS::Flavour   m_flav;
    matter_form::code m_form;
    double m_fraction1, m_radius1, m_radius2;

    void SetDefaults(const ATOOLS::Flavour & flav);
    void Initialise();
  public:
    Form_Factor(const ATOOLS::Flavour & flav);
    ~Form_Factor() {};

    ATOOLS::Vec4D operator()();
    inline const matter_form::code & GetForm() const { return m_form; }
    inline const double & Radius1()            const { return m_radius1; }
    inline const double & Radius2()            const { return m_radius2; }
    inline const double & Fraction1()          const { return m_fraction1; }
  };
}

#endif
