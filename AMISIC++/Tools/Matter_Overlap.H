#ifndef AMISIC_Tools_Matter_Overlap_H
#define AMISIC_Tools_Matter_Overlap_H

#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Vector.H"


namespace AMISIC {
  class Matter_Overlap : public ATOOLS::Function_Base {
  private:
    overlap_form::code m_overlapform;
    double m_radius1, m_radius2, m_radius3;
    double m_radius12, m_radius22, m_radius32, m_fraction1;
    double m_norm1, m_norm2, m_norm3;
    double m_bstep, m_bmax, m_norm, m_integral;
    
    void InitializeFormFactors();
    void CalculateIntegral();
 public:
    Matter_Overlap();
    ~Matter_Overlap();

    void   Initialize();    
    double operator()(double b);
    double SelectB(const bool & mode=false) const;

    ATOOLS::Vec4D SelectPositionForScatter(const double & b) const;
    ATOOLS::Vec4D SelectRelativePositionForParton() const;

    const double & Bstep()    const { return m_bstep; }
    const double & Bmax()     const { return m_bmax; }
    const double & Integral() const { return m_integral; }
    
    inline const double & Radius() const {
      switch (m_overlapform) {
      case overlap_form::Single_Gaussian:
        return m_radius1;
      case overlap_form::Double_Gaussian:
        return m_radius3;
      }
    }
  };

  class MO_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
  public:
    MO_Integrand(Matter_Overlap * mo) : p_mo(mo) {}
    ~MO_Integrand() {}
    double operator()(double b);
  };
}

#endif
