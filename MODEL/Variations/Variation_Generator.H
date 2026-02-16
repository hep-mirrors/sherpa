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
        /*
        Class to steer the generation of model parameter variations, used from Single_Process
        */
        class Variation_Generator : public Base{
            public:
                /* Constructors */
                Variation_Generator(const Args& args);                                  // constructor @param args struct containing pointer to Single_Process @returns Variation_Generator

                /* Public Functions */
                void GenerateAndFillWeightsMap(ATOOLS::Weights_Map& wgtmap) override;   // generate additional weights for an event @param wgtmap Weights_Map, this is filled with new weights @returns None
                void ResetWeightsMap(ATOOLS::Weights_Map& wgtmap) override;             // reset the weightsmap regarding these variations @param wgtmap Weights_Map, to be reset @returns None
                inline bool IsOkay() const {return okay;}                               // check if this object is valid @param None @returns bool

            protected:  
                /* Attributes */
                bool okay;                                  // validity of this object
                PHASIC::Single_Process* p_proc;             // pointer to the Single_Process this belongs to
                Variations* p_vars;                         // pointer to the variations of the model

                /* Internal Functions */
                double Calculate();                         // calculate the partonic cross section for the process this belongs to @param None @return double
                void UpdateAllCouplings(VariationKey key);  // update all model constants and couplings to the specified variation @param key VariationKey @returns None
                void SetConstants(VariationKey key);        // set all model constants to the specified variation @param key VariationKey @returns None
                void SetCouplings(VariationKey key);        // set all model couplings to the specified variation @param key VariationKey @returns None
                void ResetAllCouplings();                   // reset model constants and couplings to the nominal @param None @returns None
        };
    };
};

#endif