#include "MODEL/Variations/Variation_Generator.H"

namespace MODELVARIATIONS {
    Variation_Generator::Variation_Generator(const Args& args){
        p_proc = dynamic_cast<PHASIC::Single_Process*> (args.p_proc);
        okay = true;
        okay &= (p_proc != nullptr);
        okay &= MODEL::s_model->InitVariations();
        p_vars = MODEL::s_model->GetParameterVariations();
        okay &= (p_vars != nullptr);
        okay &= p_vars->IsOkay();
    }

    /*
    Destructor, does nothing
    */
    Variation_Generator::~Variation_Generator(){}

    /*
    generate the map for the weights under "ParameterVariations", then calculates the weights
    */
    void Variation_Generator::GenerateAndFillWeightsMap(ATOOLS::Weights_Map& wgtmap){
        double nominal = p_proc->LastXS();
        for (VariationKey var : p_vars->GetVariations()){
            if (nominal == 0.0) {
                // save some time, fix ratio
                wgtmap["ParameterVariations"][var.Identifier()] = 0;
                continue;
            }
            msg_Debugging() << "switch to " << var;
            UpdateAllCouplings(var);
            double part = Calculate();
            double weight = part/nominal;
            msg_Debugging() << " nominal: " << nominal << ", current: " << part << ", weight: " << weight << std::endl;
            wgtmap["ParameterVariations"][var.Identifier()] = weight;
        }
        ResetAllCouplings();
    }

    void Variation_Generator::UpdateAllCouplings(VariationKey key){
        // reset first
        SetConstants(p_vars->Nominal());
        SetConstants(key);
        SetCouplings(key);
    }

    void Variation_Generator::SetConstants(VariationKey key){
        msg_Debugging() << "Updating Constants to " << key.Identifier() << std::endl;
        for (size_t i = 0; i < key.Size(); i++) {
            if (MODEL::s_model->Constants()->count(key.Name(i)) == 1) MODEL::s_model->Constants()->at(key.Name(i)) = key.Value(i);
        }
    }

    void Variation_Generator::SetCouplings(VariationKey key){
        msg_Debugging() << "Updating dependent Vertices..." << std::endl;
        for (std::string name : key.Names()){
            std::set<MODEL::Single_Vertex*>* s_dependents = p_vars->Dependents(name);
            for (std::set<MODEL::Single_Vertex*>::iterator it_v = s_dependents->begin(); it_v != s_dependents->end(); it_v++){
                MODEL::Single_Vertex* p_v = *it_v;
                p_v->UpdateCouplings(MODEL::s_model->Constants());
            }
        }
    }

    void Variation_Generator::ResetAllCouplings(){
        SetConstants(p_vars->Nominal());
        SetCouplings(p_vars->Nominal());
    }

    /*
    Calculates the weight of the event under the specified variation
    */
    double Variation_Generator::Calculate(){
        return p_proc->Partonic(p_proc->Integrator()->Momenta());
    }

    /*
    Reset all weights, specifically set them to 1.0
    */
    void Variation_Generator::ResetWeightsMap(ATOOLS::Weights_Map& wgtmap){
        for (auto var : p_vars->GetVariations()){
            wgtmap["ParameterVariations"][var.Identifier()] = 1.0;
        }
    }
}

DECLARE_GETTER(MODELVARIATIONS::Variation_Generator, "MODEL_PARAMETERS", Base, Args);

Base* ATOOLS::Getter<Base, Args, MODELVARIATIONS::Variation_Generator>::
operator()(const Args& args) const
{
    MODELVARIATIONS::Variation_Generator* var_gen = new MODELVARIATIONS::Variation_Generator(args);
    if (var_gen->IsOkay()) return var_gen;
    return nullptr;
}

void ATOOLS::Getter<Base, Args, MODELVARIATIONS::Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
str << "Info for Model Param Variations \n";
}

