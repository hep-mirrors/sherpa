#ifndef PHASIC_Main_Phase_Space_Variation_Weights_H
#define PHASIC_Main_Phase_Space_Variation_Weights_H

namespace ATOOLS {
  class Variation_Weights;
}

namespace PHASIC {
  class Process_Base;
  class Phase_Space_Handler;

  class Phase_Space_Variation_Weights {
  private:
    ATOOLS::Variation_Weights * p_weights;
  public:
    inline void SetVariationWeights(ATOOLS::Variation_Weights *vw) { p_weights = vw; }
    inline ATOOLS::Variation_Weights * operator()() { return p_weights; }
  };
}

#endif
