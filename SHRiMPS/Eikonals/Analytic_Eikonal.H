#ifndef SHRIMPS_Eikonals_Analytic_Eikonal_H
#define SHRIMPS_Eikonals_Analytic_Eikonal_H

#include "SHRiMPS/Tools/Parameter_Structures.H"
#include "ATOOLS/Math/Function_Base.H"

namespace SHRIMPS {
  class Analytic_Eikonal : public ATOOLS::Function_Base {
  private:
    double m_norm, m_prefactor, m_expnorm;
  public:
    Analytic_Eikonal();
    double operator()(const double & B) const;
  };
}

#endif
