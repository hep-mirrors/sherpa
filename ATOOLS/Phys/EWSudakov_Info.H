#ifndef ATOOLS_Phys_EWSudakov_Info_H
#define ATOOLS_Phys_EWSudakov_Info_H

#include <iostream>

namespace ATOOLS {

  enum class EWSudakov_Log_Type {
    // TODO: add more coeffs (and their implementation)
    Ls,       ///< correction \propto \log(s/m_mW2)^2
    lZ,       ///< correction \propto \log(s/m_mW2)
    lSSC,     ///< correction \propto \log(s/m_mW2) \log(|r_kl|/s)
    lC,       ///< collinear non-Yukawa corrections \propto \log(s/m_mW2)
    lYuk,     ///< collinear Yukawa corrections \propto \log(s/m_mW2)
    lPR,      ///< parameter renormalisation corrections \propto \log(s/m_mW2)
  };

  EWSudakov_Log_Type EWSudakovLogTypeFromString(const std::string&);
  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type&);

}

#endif

