#ifndef THRew__Variation_Generator__H
#define THRew__Variation_Generator__H

#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "AddOns/THRew/KFactor.H"

namespace THRew {

  class Variation_Generator
      : public ATOOLS::Hard_Process_Variation_Generator_Base {
  public:

    Variation_Generator(
        const ATOOLS::Hard_Process_Variation_Generator_Arguments&);

    void GenerateAndFillWeightsMap(ATOOLS::Weights_Map&) override;
    void ResetWeightsMap(ATOOLS::Weights_Map&) override;

  private:

    THRew_KFactor m_kfactor;
  };

} // namespace THRew

#endif
