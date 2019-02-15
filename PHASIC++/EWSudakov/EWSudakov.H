#ifndef PHASIC_EWSudakov_EWSudakov_H
#define PHASIC_EWSudakov_EWSudakov_H

#include "ATOOLS/Phys/Flavour.H"

#include <array>
#include <set>
#include <unordered_map>

/// definitions for types and enums used in the module
/** Note that this module refers in many places to a publication by
 *  Denner/Pozzorini, which can be found under arXiv:hep-ph/0010201.
 *  Also the definitions and calculations follow their conventions
 *  as far as possible to make comparisons and checks more straightforward.
 */

namespace PHASIC {

  using Two_Leg_Indizes = std::array<size_t, 2>;
  using Leg_Index_Set = std::set<size_t>;
  using Leg_Set = std::set<std::pair<size_t, kf_code>>;
  using LegIndex_Flavour_Map = std::unordered_map<size_t, ATOOLS::Flavour>;
  using Couplings = std::unordered_map<kf_code, double>;

  enum class EWSudakov_Log_Type {
    // TODO: add more coeffs (and their implementation)
    Ls,       ///< correction \propto \log(s/m_mW2)^2
    lZ,       ///< correction \propto \log(s/m_mW2)
    lSSC,     ///< correction \propto \log(s/m_mW2) \log(|r_kl|/s)
    lC,       ///< collinear non-Yukawa corrections \propto \log(s/m_mW2)
    lYuk,     ///< collinear Yukawa corrections \propto \log(s/m_mW2)
    lPR,      ///< parameter renormalisation corrections \propto \log(s/m_mW2)
  };

  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type&);
  std::ostream& operator<<(std::ostream& os, const Leg_Set&);

}

#endif
