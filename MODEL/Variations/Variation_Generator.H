#ifndef MODEL_Variations_Variation_Generator_H
#define MODEL_Variations_Variation_Generator_H

#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"

using Base = ATOOLS::Hard_Process_Variation_Generator_Base;
using Args = ATOOLS::Hard_Process_Variation_Generator_Arguments;

namespace MODEL{
    namespace VARIATIONS {
        class Variation_Generator : public Base{
            public:
                Variation_Generator(const Args& args);
                ~Variation_Generator() override;
                void GenerateAndFillWeightsMap(ATOOLS::Weights_Map&) override;
                void ResetWeightsMap(ATOOLS::Weights_Map&) override;
                inline bool IsOkay() const {return okay;}

            protected:  
                bool okay;
                PHASIC::Single_Process* p_proc;
                Variations* p_vars;

                double Calculate();            
                void UpdateAllCouplings(VariationKey key);
                void SetConstants(VariationKey key);
                void SetCouplings(VariationKey key);
                
                void ResetAllCouplings();
        };
    };
};

#endif