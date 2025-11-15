#ifndef AddOns_UFOVariations_Variation_Generator_H
#define AddOns_UFOVariations_Variation_Generator_H

#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "AddOns/UFOVariations/UFO_Variation_Key.H"
#include "AddOns/UFOVariations/Variations.H"
#include <vector>
#include <map>

using Base = ATOOLS::Hard_Process_Variation_Generator_Base;
using Args = ATOOLS::Hard_Process_Variation_Generator_Arguments;

namespace UFOVariations{
    class Variation_Generator : public Base{
        public:
            Variation_Generator(const Args& args);
            ~Variation_Generator() override;
            void GenerateAndFillWeightsMap(ATOOLS::Weights_Map&) override;
            void ResetWeightsMap(ATOOLS::Weights_Map&) override;
            void UpdateAllCouplings(VariationKey key);

        protected:
            PHASIC::Single_Process* p_proc;
            Variations* p_vars;
            double Calculate();
    };
};

#endif