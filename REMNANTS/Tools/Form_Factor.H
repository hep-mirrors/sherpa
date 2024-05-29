#ifndef REMNANTS_FORM_FACTOR_H
#define REMNANTS_FORM_FACTOR_H

#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"

namespace REMNANTS {
  class Form_Factor {
  private:
    ATOOLS::Flavour   m_flav;
    matter_form       m_form;
    double m_fraction1, m_radius1, m_radius2, m_softexp;

    void SetDefaults(const ATOOLS::Flavour & flav);
    void Initialise();
  public:
    Form_Factor(const ATOOLS::Flavour & flav);
    ~Form_Factor() {};

    ATOOLS::Vec4D operator()(const double & x=-1., const double & Q2=-1.);
    double        B(const double & x=-1., const double & Q2=-1.);
    const double  Radius(const double & x=-1.,const double & Q2=-1.) const;
    inline const matter_form & GetForm()      const { return m_form; }
    inline const double      & Radius1()      const { return m_radius1; }
    inline const double      & Radius2()      const { return m_radius2; }
    inline const double      & Fraction1()    const { return m_fraction1; }
    inline const double      & SoftExponent() const { return m_softexp; }
  };
}

#endif
