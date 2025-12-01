#include "MODEL/UFO/Variation_Generator.H"

namespace UFO {
    Variation_Generator::Variation_Generator(const Args& args){
        p_proc = dynamic_cast<PHASIC::Single_Process*> (args.p_proc);
        if (!p_proc) THROW(fatal_error, "No Single Process was given for Variation.")
        UFO_Model* model = dynamic_cast<UFO_Model*> (MODEL::s_model);
        if (!model) THROW(fatal_error, "The model does not seem to implement Variations of external parameters :(")
        p_vars = model->GetParameterVariations();
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
            msg_Debugging() << "switch to " << var;
            UpdateAllCouplings(var);
            double part = Calculate();
            msg_Debugging() << " nominal: " << nominal << ", current: " << part << std::endl;
            wgtmap["UFOVariations"][var.Identifier()] = part/nominal;
        }
        // reset to default vertices TODO save nominal param values somewhere
        UpdateAllCouplings(p_vars->Nominal());
        p_proc->SetLookUp(true);
    }

    void Variation_Generator::UpdateAllCouplings(VariationKey key){
            // reset p_consts
            msg_Debugging() << "Resetting Constants..." << std::endl;
            VariationKey nominal = p_vars->Nominal();
            for (size_t i = 0; i < nominal.Size(); i++) {
                if (MODEL::s_model->Constants()->count(nominal.Name(i)) == 1) MODEL::s_model->Constants()->at(nominal.Name(i)) = nominal.Value(i);
            }
            // update p_consts
            msg_Debugging() << "Updating Constants to " << key.Identifier() << std::endl;
            for (size_t i = 0; i < key.Size(); i++) {
                if (MODEL::s_model->Constants()->count(key.Name(i)) == 1) MODEL::s_model->Constants()->at(key.Name(i)) = key.Value(i);
            }
            // update Kabbalas
            msg_Debugging() << "Updating dependent Vertices..." << std::endl;
            for (std::string name : key.Names()){
                std::set<MODEL::Single_Vertex*>* s_dependents = p_vars->Dependents(name);
                for (std::set<MODEL::Single_Vertex*>::iterator it_v = s_dependents->begin(); it_v != s_dependents->end(); it_v++){
                    MODEL::Single_Vertex* p_v = *it_v;
                    p_v->UpdateCouplings(MODEL::s_model->Constants());
                }
                   
            }
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
            wgtmap["UFOVariations"][var.Identifier()] = 1.0;
        }
    }
}

DECLARE_GETTER(UFO::Variation_Generator, "UFOVariations", Base, Args);

Base* ATOOLS::Getter<Base, Args, UFO::Variation_Generator>::
operator()(const Args& args) const
{
    return new UFO::Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, UFO::Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
str << "Info for UFO Param Variations \n";
}

