#ifndef PHASIC_EWSudakov_EWSudakov_H
#define PHASIC_EWSudakov_EWSudakov_H

#include <array>
#include <set>

/// definitions for types and enums used in the module

namespace PHASIC {

  using Two_Leg_Indizes = std::array<size_t, 2>;
  using Leg_Index_Set = std::set<size_t>;

  enum class EWSudakov_Log_Type {
    // TODO: add more coeffs (and their implementation)
    Ls,       ///< correction \propto \log(s/m_mW2)^2
    lZ,       ///< correction \propto \log(s/m_mZ2)
    lSSC,     ///< correction \propto \log(s/m_mW2) \log(|r_kl|/s)
  };

  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type&);

}

#endif
