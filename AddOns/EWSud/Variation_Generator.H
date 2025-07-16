#ifndef AddOns_EWSud_Variation_Generator_H
#define AddOns_EWSud_Variation_Generator_H

#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "AddOns/EWSud/KFactor.H"

namespace EWSud {

  class Variation_Generator :
    public ATOOLS::Hard_Process_Variation_Generator_Base {

  public:

    Variation_Generator(
        const ATOOLS::Hard_Process_Variation_Generator_Arguments&);

    void GenerateAndFillWeightsMap(ATOOLS::Weights_Map&) override;

    void ResetWeightsMap(ATOOLS::Weights_Map&) override;

  private:

    Sudakov_KFactor m_kfactor;
    bool m_ewsudakov_rs;

  };

}

#endif
