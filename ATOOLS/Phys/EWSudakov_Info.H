#ifndef ATOOLS_Phys_EWSudakov_Info_H
#define ATOOLS_Phys_EWSudakov_Info_H

#include <iostream>
#include <map>

namespace ATOOLS {

  enum class EWSudakov_Log_Type {
    Ls,       ///< correction \propto \log(s/m_mW2)^2
    lZ,       ///< correction \propto \log(s/m_mW2)
    lSSC,     ///< correction \propto \log(s/m_mW2) \log(|r_kl|/s)
    lC,       ///< collinear non-Yukawa corrections \propto \log(s/m_mW2)
    lYuk,     ///< collinear Yukawa corrections \propto \log(s/m_mW2)
    lPR,      ///< parameter renormalisation corrections \propto \log(s/m_mW2)
  };

  EWSudakov_Log_Type EWSudakovLogTypeFromString(const std::string&);
  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type&);

  class EWSudakov_Log_Corrections_Map
      : public std::map<EWSudakov_Log_Type, double> {
  public:
    double KFactor() const;
  };

  struct Mandelstam_Variables {
    Mandelstam_Variables(double _s, double _t, double _u):
      s {_s}, t {_t}, u {_u}
    {};
    double s {0.0};
    double t {0.0};
    double u {0.0};
  };

}

#endif
