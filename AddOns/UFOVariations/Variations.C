#include "AddOns/UFOVariations/Variations.H"

#include "AddOns/UFOVariations/Combinations.H"
#include "ATOOLS/Org/Settings.H"
#include "PHASIC++/Process/Single_Process.H"
#include <ostream>
#include <map>
#include <algorithm>

#define MAX_VARIATION_NUMBER 1000

namespace UFOVariations {
    /*
    Constructs the Variations Object
    Iterates over the "UFO_VARIATIONS" Setting and adds the Variations to its own list
    */
    Variations::Variations(){
        // Deal with UFO Param Variations
        // get settings and read them to map
        ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
        std::vector<ATOOLS::Scoped_Settings> items = s["UFO_VARIATIONS"].GetItems();
        std::map<std::string, std::vector<double_t>> variation_map = std::map<std::string, std::vector<double_t>>();
        for (ATOOLS::Scoped_Settings item : items) ReadSingleParamVariation(item, &variation_map);
        // map is done here, add variation key elements from combinations
        Combinations *comb = new Combinations(&variation_map); 
        comb->AddAllCombinations(&v_variations);
        // check if there are too many variations
        if (Size() >= MAX_VARIATION_NUMBER) THROW(normal_exit, "You are trying too many Variations, please reconsider.");
        // Register the dependent Kabbalas of the Model vertices
        FindDependentCouplings();
        StoreNominal();
    }

    /*
    Destructor, does nothing
    */
    Variations::~Variations(){}

    /*
    read a single param variation (one variable, mulitple values) and put it in a map
    */
    void Variations::ReadSingleParamVariation(ATOOLS::Scoped_Settings& s, std::map<std::string, std::vector<double_t>> *variation_map){
        // read and check variable name
        std::string name = s["NAME"].SetDefault("None").Get<std::string>();
        variables.push_back(name);
        if (name == "None") THROW(invalid_input, "Variable name must be given for UFO Variation");
        // read in values
        std::vector<double_t> values = s["VALUES"].SetDefault("0").GetVector<double_t>();
        // TODO :: check values
        // save to map
        variation_map->insert(std::make_pair(name, values));
    }

    // singleton init
    void Variations::Init(){
        if(variations == nullptr){
            msg_Out() << "Looking for UFO Param Variations..." << std::endl;
            variations = new Variations();
        }
    }

    // singleton getter
    Variations* Variations::Get() {
        Variations::Init();
        return variations;
    }

    // find and register the dependencies of the couplings of all Model vertices
    void Variations::FindDependentCouplings() {
        // for all the varied parameters
        for (std::string var_name : variables){
            // empty vector init
            dependent_couplings.insert(std::make_pair(var_name, std::vector<ATOOLS::Kabbala*>()));
            // go through vertices
            for (MODEL::Single_Vertex v : MODEL::s_model->Vertices()) {
                // go through couplings of the vertex and store a pointer to it if its dependent
                for (ATOOLS::Kabbala k : v.cpl) {
                    if (k.DependsOn(var_name)) dependent_couplings[var_name].push_back(&k);
                }
            }
        }
    }

    void Variations::StoreNominal() {
        std::vector<double_t> values = {};
        for (std::string name : variables){
            values.push_back(MODEL::s_model->Constants()->at(name));
        }
        nominal = VariationKey(variables, values);
    }

    // singleton
    Variations* Variations::variations = nullptr;
}

