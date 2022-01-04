#ifndef EWSudakov_EWSudakov_H
#define EWSudakov_EWSudakov_H

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

namespace EWSud {

  enum class EWSudakov_Log_Type {
    Ls,       ///< correction \propto \log(s/m_mW2)^2
    lZ,       ///< correction \propto \log(s/m_mW2)
    lSSC,     ///< correction \propto \log(s/m_mW2) \log(|r_kl|/s)
    lC,       ///< collinear non-Yukawa corrections \propto \log(s/m_mW2)
    lYuk,     ///< collinear Yukawa corrections \propto \log(s/m_mW2)
    lPR,      ///< parameter renormalisation corrections \propto \log(s/m_mW2)
    lI,       ///< I operator corresponding to clustered pairs
  };

  EWSudakov_Log_Type EWSudakovLogTypeFromString(const std::string&);
  std::ostream& operator<<(std::ostream& os, const EWSudakov_Log_Type&);

  class EWSudakov_Log_Corrections_Map
      : public std::map<EWSudakov_Log_Type, double> {
  public:
    double KFactor() const;
    friend std::ostream& operator<<(std::ostream&,
                                    const EWSudakov_Log_Corrections_Map&);
  };

  struct Mandelstam_Variables {
    Mandelstam_Variables(double _s, double _t, double _u):
      s {_s}, t {_t}, u {_u}
    {};
    double s {0.0};
    double t {0.0};
    double u {0.0};
  };

  /** High energy schemes:
  *
  * _default: DP-scheme, that is do not compute the kfactor if at least one of the
  *          invariants isn't large enough
  *
  * _tolerant: Only discard those pair of particles whose invariant isn't large
  *            enough from the calculation
  *
  * _cluster_dumb: if a pair of particles has an insufficient invariant mass,
  *               try combining the two momenta and with the combined momentum
  *               do the same as default
  *
  * _cluster_l1: cluster particles that can form a resonant decay
  *
  */
  enum class EWSudakovHighEnergySchemes{
    _default,
    _tolerant,
    _cluster_dumb,
    _cluster_l1,
  };

  using Two_Leg_Indizes = std::array<size_t, 2>;
  using Leg_Index_Set = std::set<size_t>;
  // NOTE: the following map can not be unordered, since it is being used as a
  // key and therefore needs to support comparison operations
  using Leg_Kfcode_Map = std::map<size_t, kf_code>;
  // NOTE: Flavour has a (long int)-conversion operator
  using Leg_Kfcode_Map_Signed = std::unordered_map<size_t, long int>;
  using Couplings = std::unordered_map<long int, Complex>;

  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map&);
  std::ostream& operator<<(std::ostream& os, const Leg_Kfcode_Map_Signed&);

  Leg_Kfcode_Map ConvertToPhysicalPhase(Leg_Kfcode_Map legs);

}

#endif
