#ifndef AMISIC_Tools_Matter_Overlap_H
#define AMISIC_Tools_Matter_Overlap_H

#include "REMNANTS/Main/Remnant_Handler.H"
#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Vector.H"


namespace AMISIC {
  class Matter_Overlap : public ATOOLS::Function_Base {
  private:
    overlap_form m_overlapform;
    double m_bstep, m_bmax, m_norm, m_integral;
    double m_radius[4], m_radius2[4], m_rnorm[4], m_fraction[4];

    void InitializeFormFactors(REMNANTS::Remnant_Handler * remnant_handler);
    void CalculateIntegral();
 public:
    Matter_Overlap();
    ~Matter_Overlap();

    void   Initialize(REMNANTS::Remnant_Handler * remnant_handler);
    double operator()(double b);
    double SelectB(const bool & mode=false) const;

    ATOOLS::Vec4D SelectPositionForScatter(const double & b) const;
    ATOOLS::Vec4D SelectRelativePositionForParton() const;

    const double & Bmax()     const { return m_bmax; }
    const double & Integral() const { return m_integral; }
  };

  class MO_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
  public:
    MO_Integrand(Matter_Overlap * mo) : p_mo(mo) {}
    ~MO_Integrand() = default;
    double operator()(double b);
  };
}

#endif
