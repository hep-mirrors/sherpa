#ifndef METOOLS_HadronCurrents_A0_Decays_H
#define METOOLS_HadronCurrents_A0_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  // FIXME: locally-defined kf-code constants for a0(980)/a0(1450),
  // since this tarball's copy of ATOOLS/Phys/Flavour_Tags.H doesn't
  // have them yet - see the identical note in Kstar_Decays.H/K1_Decays.H.
  // a0(980)/a0(1450): kf-codes confirmed against the real
  // ATOOLS/Phys/Flavour_Tags.H. Note the REAL macro names use an
  // underscore ("a_0", not "a0") - kf_a_0_980(_plus)=9000111/9000211,
  // kf_a_0_1450(_plus)=10111/10211 (values match what was guessed
  // here, but the earlier "kf_a0_..." names did not exist and would
  // not have compiled). Use the real names directly - no local
  // aliases declared, since they are #define macros in
  // Flavour_Tags.H (a local `const kf_code kf_a_0_980 = ...;` would
  // collide with the macro substitution).

  // a0(980), scalar isotriplet resonance (kf 9000111/9000211), decaying
  // to eta pi (100% when eta' pi is closed) plus eta' pi with the SAME
  // coupling once that channel opens (Frank, private communication) -
  // see A0_Decays.C for how "same coupling" is implemented technically.
  class A0_980_0_Lineshape : public Total_Width_Base {
  public:
    A0_980_0_Lineshape();
  };

  class A0_980_plus_Lineshape : public Total_Width_Base {
  public:
    A0_980_plus_Lineshape();
  };

  // a0(1450), scalar isotriplet resonance (kf 10111/10211), decaying
  // into {pi eta = 10%, pi eta' = 3%, K Kbar = 10%, omega pi pi = 67%}
  // (Frank, private communication). The omega pi pi channel has no
  // dedicated 3-body phase-space treatment in this code - see the
  // Constant_PW flag in A0_Decays.C.
  class A0_1450_0_Lineshape : public Total_Width_Base {
  public:
    A0_1450_0_Lineshape();
  };

  class A0_1450_plus_Lineshape : public Total_Width_Base {
  public:
    A0_1450_plus_Lineshape();
  };
}

#endif
