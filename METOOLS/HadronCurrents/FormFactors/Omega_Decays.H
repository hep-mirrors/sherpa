#ifndef METOOLS_HadronCurrents_Omega_Decays_H
#define METOOLS_HadronCurrents_Omega_Decays_H

#include "METOOLS/HadronCurrents/FormFactors/Resonance_Base.H"

namespace METOOLS {
  class Omega_782_Lineshape : public Total_Width_Base {
  public:
    Omega_782_Lineshape();
  };
  class Omega_1420_Lineshape : public Total_Width_Base {
  public:
    Omega_1420_Lineshape();
  };
  class Omega_1600_Lineshape : public Total_Width_Base {
  public:
    Omega_1600_Lineshape();
  };

  // phi(1020), kf 333, decaying 50/50 into K+K- and K0 K0bar (Frank,
  // private communication).
  class Phi_1020_Lineshape : public Total_Width_Base {
  public:
    Phi_1020_Lineshape();
  };
}

#endif
