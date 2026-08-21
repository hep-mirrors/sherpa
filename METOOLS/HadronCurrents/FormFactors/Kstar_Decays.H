#ifndef METOOLS_HadronCurrents_Kstar_Decays_H
#define METOOLS_HadronCurrents_Kstar_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  // FIXME: locally-defined kf-code constants for K*_0(1430), since
  // this tarball's copy of ATOOLS/Phys/Flavour_Tags.H doesn't have
  // them yet. Please add proper named aliases there if you'd like
  // other code to reference these particles by name - defining them
  // here (in a header, so all translation units see the SAME values)
  // is enough for the numerics in this patch to work regardless.
  // K*_0(1430) and K*(1680): kf-codes now confirmed against the real
  // ATOOLS/Phys/Flavour_Tags.H (kf_K_0_star_1430(_plus)=10311/10321,
  // kf_K_star_1680(_plus)=30313/30323) - use those directly, no local
  // aliases needed (and none should be declared here, since
  // Flavour_Tags.H defines them as #define macros: a local `const
  // kf_code kf_K_star_1680 = ...;` would collide with the macro
  // substitution and fail to compile once that header is in scope).

  class Kstar_892_0_Lineshape : public Total_Width_Base {
  public:
    Kstar_892_0_Lineshape();
  };

  class Kstar_892_plus_Lineshape : public Total_Width_Base {
  public:
    Kstar_892_plus_Lineshape();
  };

  class Kstar_1410_0_Lineshape : public Total_Width_Base {
  public:
    Kstar_1410_0_Lineshape();
  };

  class Kstar_1410_plus_Lineshape : public Total_Width_Base {
  public:
    Kstar_1410_plus_Lineshape();
  };

  // K*_0(1430), scalar kaon resonance (kf 10311/10321), decaying into
  // K pi only - see Kstar_Decays.C for the isospin split.
  class Kstar0_1430_0_Lineshape : public Total_Width_Base {
  public:
    Kstar0_1430_0_Lineshape();
  };

  class Kstar0_1430_plus_Lineshape : public Total_Width_Base {
  public:
    Kstar0_1430_plus_Lineshape();
  };

  // K*(1680), vector kaon resonance (kf 30313/30323), decaying in
  // equal (1/3, 1/3, 1/3) parts to K pi, K rho(770), K*(892) pi
  // (Frank, private communication). rho(770) and K*(892) are treated
  // as on-shell in the K rho / K* pi channels (V_VP, same treatment as
  // K1(1270)/K1(1400) in K1_Decays.C), not integrated over their own
  // lineshapes.
  class Kstar_1680_0_Lineshape : public Total_Width_Base {
  public:
    Kstar_1680_0_Lineshape();
  };

  class Kstar_1680_plus_Lineshape : public Total_Width_Base {
  public:
    Kstar_1680_plus_Lineshape();
  };
}

#endif
