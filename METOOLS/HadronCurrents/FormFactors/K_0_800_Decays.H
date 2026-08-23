#ifndef METOOLS_HadronCurrents_K_0_800_Decays_H
#define METOOLS_HadronCurrents_K_0_800_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  // K_0*(700)/kappa - the light strange scalar, completing the light
  // scalar nonet alongside f0(600) (F0_Decays.H) and a0(980)
  // (A0_Decays.H). kf-codes 9000311 (neutral)/9000321 (charged) -
  // NEWLY added (not previously in ATOOLS/Phys/Flavour_Tags.H), same
  // "9000XXX = light nonet" numbering pattern already used for
  // f0(600)=9000221 and a0(980)=9000111/9000211, paralleling how
  // K_0*(1430) (heavier nonet partner) sits at 10311/10321 alongside
  // f0(1370)=10221 and a0(1450)=10111/10211.
  //
  // Decay content: essentially 100% K pi (kappa sits well below any
  // other open channel), split 50/50 by isospin - identical structure
  // to Kstar0_1430_0/plus_Lineshape (Kstar_Decays.H/.C), just at a
  // different mass/width. Uses the same S_PP scalar partial width.
  //
  // CAVEAT: kappa's pole position is notoriously poorly determined
  // experimentally - quoted values range roughly M~630-900 MeV,
  // Gamma~250-650 MeV depending on the analysis/extraction method,
  // even more unsettled than f0(600)'s own (also broad, also
  // model-dependent) parameters. The mass/width used here should be
  // treated as an explicit, overridable placeholder - see
  // K_0_800_Decays.C for the specific numbers and their (lack of a
  // single authoritative) source.
  class K_0_800_0_Lineshape : public Total_Width_Base {
  public:
    K_0_800_0_Lineshape();
  };

  class K_0_800_plus_Lineshape : public Total_Width_Base {
  public:
    K_0_800_plus_Lineshape();
  };
}

#endif
