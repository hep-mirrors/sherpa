#ifndef SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#define SHRIMPS_Cross_Sections_Sigma_Inelastic_H
#include "SHRiMPS/Cross_Sections/Sigma_Base.H"
#include <map>

namespace SHRIMPS {
  class Sigma_Inelastic : public Sigma_Base {
  private:
  public:
    Sigma_Inelastic() {}
    ~Sigma_Inelastic() {}

    double GetValue(const double & B); 
    double GetCombinedValue(const double & B);
    double XSec() { return m_sigma; }

    std::vector<double> * FillBGrid(Omega_ik * eikonal);    

    double Test();
 };
}
#endif
