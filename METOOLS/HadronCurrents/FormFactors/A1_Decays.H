#ifndef METOOLS_HadronCurrents_A1_Decays_H
#define METOOLS_HadronCurrents_A1_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  class A1_1260_0_Lineshape : public Total_Width_Base {
  public:
    A1_1260_0_Lineshape();
  };

  class A1_1260_plus_Lineshape : public Total_Width_Base {
  public:
    A1_1260_plus_Lineshape();
  };

  // Running width for the a1(1260)+ as used in the KS form factor:
  // pure rho(770) only in the Dalitz integral, no rho(1450) or f0 admixture.
  class KS_A1_1260_plus_Width : public Total_Width_Base {
  public:
    KS_A1_1260_plus_Width();
  };
}

#endif
