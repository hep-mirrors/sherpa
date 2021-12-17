#ifndef AddOns_EWSud_Variation_Weight_Names_H
#define AddOns_EWSud_Variation_Weight_Names_H

#include "ATOOLS/Phys/Hard_Process_Variation_Weight_Names.H"
#include "AddOns/EWSud/KFactor.H"

namespace EWSud {

  class Variation_Weight_Names :
    public ATOOLS::Hard_Process_Variation_Weight_Names_Base {

  public:

    Variation_Weight_Names(
        const ATOOLS::Hard_Process_Variation_Weight_Names_Arguments&);

    std::string WeightMapKey() const override;
    std::string WeightNameForVariation(const std::string&) const override;

  };

}

#endif
