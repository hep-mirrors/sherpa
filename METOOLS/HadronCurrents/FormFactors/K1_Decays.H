#ifndef METOOLS_HadronCurrents_K1_Decays_H
#define METOOLS_HadronCurrents_K1_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  // FIXME: locally-defined kf-code constants for K1(1270)/K1(1400),
  // since this tarball's copy of ATOOLS/Phys/Flavour_Tags.H doesn't
  // have them yet - see the identical note in Kstar_Decays.H for
  // K*_0(1430). Defined here (in a header) so every translation unit
  // that includes this file (including FF_0_PPP.C) sees the same
  // values.
  // K1(1270)/K1(1400): kf-codes confirmed against the real
  // ATOOLS/Phys/Flavour_Tags.H (kf_K_1_1270(_plus)=10313/10323,
  // kf_K_1_1400(_plus)=20313/20323, exactly matching what was guessed
  // here) - use those directly. No local aliases declared (they are
  // #define macros in Flavour_Tags.H; a local `const kf_code
  // kf_K_1_1270 = ...;` would collide with the macro substitution).

  // K1(1270), axial-vector kaon resonance (kf 10313/10323), decaying
  // into on-shell rho(770)K and K*(892)pi at a 70:30 mixture (Frank,
  // private communication) - see K1_Decays.C for the isospin split.
  class K1_1270_0_Lineshape : public Total_Width_Base {
  public:
    K1_1270_0_Lineshape();
  };

  class K1_1270_plus_Lineshape : public Total_Width_Base {
  public:
    K1_1270_plus_Lineshape();
  };

  // K1(1400), axial-vector kaon resonance (kf 20313/20323), decaying
  // exclusively into K*(892)pi (Frank, private communication).
  class K1_1400_0_Lineshape : public Total_Width_Base {
  public:
    K1_1400_0_Lineshape();
  };

  class K1_1400_plus_Lineshape : public Total_Width_Base {
  public:
    K1_1400_plus_Lineshape();
  };

  // Flatte line shape for the K1(1270), for form-factor model 106.  Standalone:
  // it takes no Total_Width_Base and reads only stable-hadron masses, so it is
  // independent of the K1 width registry above.  See K1_Decays.C.
  Propagator_Base * MakeK1_1270_Flatte();
}

#endif
