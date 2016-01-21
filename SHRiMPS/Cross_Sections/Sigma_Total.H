#ifndef SHRIMPS_Cross_Sections_Sigma_Total_H
#define SHRIMPS_Cross_Sections_Sigma_Total_H
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"

namespace SHRIMPS {
  class Sigma_Tot : public Sigma_Base {
  public:
    Sigma_Tot() {}
    double GetValue(const double & B);
    double GetCombinedValue(const double & B);
    double Test();
  };

  class Elastic_Slope : public Sigma_Base {
    double m_stot;
  public:
    Elastic_Slope(const double & stot) : m_stot(stot) {}
    double GetValue(const double & B);
    double GetCombinedValue(const double & B);
  };
}
#endif
