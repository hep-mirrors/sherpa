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

}

#endif
