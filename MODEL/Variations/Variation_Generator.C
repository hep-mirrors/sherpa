#include "MODEL/Variations/Variation_Generator.H"

namespace MODEL {
    namespace VARIATIONS {
        Variation_Generator::Variation_Generator(const Args& args){
            p_proc = dynamic_cast<PHASIC::Single_Process*> (args.p_proc);
            okay = true;
            okay &= (p_proc != nullptr);
            okay &= MODEL::s_model->InitVariations();
            p_vars = MODEL::s_model->GetParameterVariations();
            okay &= (p_vars != nullptr);
            if (okay) okay &= p_vars->IsOkay();
        }

        void Variation_Generator::GenerateAndFillWeightsMap(ATOOLS::Weights_Map& wgtmap){
            // nominal cs to compare
            double nominal = p_proc->LastXS();
            for (VariationKey var : p_vars->GetVariations()){
                if (nominal == 0.0) {
                    // save some time, fix ratio
                    wgtmap["ParameterVariations"][var.Identifier()] = 0;
                    continue;
                }
                // update couplings
                msg_Debugging() << "switch to " << var;
                UpdateAllCouplings(var);
                // recalculate cross section
                double part = Calculate();
                // weight is the ratio of new cs and nominal cs
                double weight = part/nominal;
                msg_Debugging() << " nominal: " << nominal << ", current: " << part << ", weight: " << weight << std::endl;
                wgtmap["ParameterVariations"][var.Identifier()] = weight;
            }
            // reset for next nominal
            ResetAllCouplings();
        }

        void Variation_Generator::UpdateAllCouplings(VariationKey key){
            // reset first
            ResetAllCouplings();
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

        double Variation_Generator::Calculate(){
            return p_proc->Partonic(p_proc->Integrator()->Momenta());
        }

        void Variation_Generator::ResetWeightsMap(ATOOLS::Weights_Map& wgtmap){
            for (auto var : p_vars->GetVariations()){
                wgtmap["ParameterVariations"][var.Identifier()] = 1.0;
            }
        }
    }
}

DECLARE_GETTER(MODEL::VARIATIONS::Variation_Generator, "MODEL_PARAMETERS", Base, Args);

/*
Get this object, but only if it is okay, return nullptr otherwise
*/
Base* ATOOLS::Getter<Base, Args, MODEL::VARIATIONS::Variation_Generator>::operator()(const Args& args) const
{
    MODEL::VARIATIONS::Variation_Generator* var_gen = new MODEL::VARIATIONS::Variation_Generator(args);
    if (var_gen->IsOkay()) return var_gen;
    return nullptr;
}

/*
Print Information on this, TODO
*/
void ATOOLS::Getter<Base, Args, MODEL::VARIATIONS::Variation_Generator>::PrintInfo(std::ostream& str, const size_t width) const
{ 
str << "Info for Model Param Variations TODO \n";
}

