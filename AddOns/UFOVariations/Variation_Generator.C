#include "AddOns/UFOVariations/Variation_Generator.H"

#include "AddOns/UFOVariations/Combinations.H"
#include "ATOOLS/Org/Settings.H"
#include "PHASIC++/Process/Single_Process.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Kabbala.H"
#include <ostream>
#include <map>
#include <algorithm>


namespace UFOVariations {
    Variation_Generator::Variation_Generator(const Args& args){
        if (!args.p_single) THROW(invalid_input, "Something is wrong, no single process");
        p_proc = args.p_single;
        p_vars = Variations::Get();
    }

    /*
    Destructor, does nothing
    */
    Variation_Generator::~Variation_Generator(){}

    /*
    generate the map for the weights under "UFOVariations", then calculates the weights
    */
    void Variation_Generator::GenerateAndFillWeightsMap(ATOOLS::Weights_Map& wgtmap){
        double nominal = p_proc->LastXS();
        p_proc->SetLookUp(false);
        for (VariationKey var : p_vars->GetVariations()){
            UpdateAllCouplings(var);
            double part = Calculate();
            msg_Out() << "swichted to " << var.Identifier() << ", nominal: " << nominal << ", current: " << part << std::endl;
            wgtmap["UFOVariations"][var.Identifier()] = part/nominal;
        }
        // reset to default vertices TODO save nominal param values somewhere
        UpdateAllCouplings(p_vars->Nominal());
        p_proc->SetLookUp(true);
    }

    void Variation_Generator::UpdateAllCouplings(VariationKey key){
            // update p_consts
            for (size_t i = 0; i < key.Size(); i++) {
                MODEL::s_model->Constants()->at(key.Name(i)) = key.Value(i);
            }
            // update Kabbalas
            for (std::string name : key.Names()){
                for (ATOOLS::Kabbala* k : p_vars->Dependents(name))
                    k->Update(MODEL::s_model->Constants());
            }
    }

    /*
    Calculates the weight of the event under the specified variation
    */
    double Variation_Generator::Calculate(){
        // TODO figure out how to redo calculation here
        ATOOLS::Vec4D_Vector vec = p_proc->Last4DVectors();
        if (vec.empty()) {
            msg_Out() << "Empty Vector for recalculating ME, returning weight 1" << std::endl;
            return 1;
        }
        return p_proc->Partonic(vec);
    }

    /*
    Reset all weights, specifically set them to 1.0
    */
    void Variation_Generator::ResetWeightsMap(ATOOLS::Weights_Map& wgtmap){
        for (auto var : p_vars->GetVariations()){
            wgtmap["UFOVariations"][var.Identifier()] = 1.0;
        }
    }
}

DECLARE_GETTER(UFOVariations::Variation_Generator, "UFOVariations", Base, Args);

Base* ATOOLS::Getter<Base, Args, UFOVariations::Variation_Generator>::
operator()(const Args& args) const
{
    return new UFOVariations::Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, UFOVariations::Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
str << "Info for UFO Param Variations \n";
}

