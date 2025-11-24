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
        msg_Debugging() << "Variable Names of Variations are: " << variables << std::endl;
        // map is done here, add variation key elements from combinations
        Combinations *comb = new Combinations(&variation_map); 
        comb->AddAllCombinations(&v_variations);
        if (msg_LevelIsDebugging()) {
            msg_Out() << "All read Variations are: " << std::endl;
            for (VariationKey key : v_variations){
                msg_Out() << "\t" << key << std::endl;
            }
        }
        else {
            msg_Out() << "\t--> Found " << Size() << " variations" << std::endl;
        }
        // check if there are too many variations
        if (Size() >= MAX_VARIATION_NUMBER) THROW(invalid_input, "You are trying too many Variations, please reconsider.");
        // Register the dependent Kabbalas of the Model vertices
        FindDependentVertices();
        StoreNominal();
        msg_Out() << "UFO Variations read." << std::endl << std::endl;
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
        std::vector<double_t> values = s["VALUES"].SetDefault("-1").GetVector<double_t>();
        if (values.back() == -1.) {
            std::vector<double_t> range = s["RANGE"].SetDefault("-1").GetVector<double_t>();
            if (range.back() == -1.) THROW(invalid_input, "Please specify concrete parameter values with VALUES: ... or a range with RANGE: [from, to, nr]");
            if (range.size()!=3) THROW(invalid_input, "Please specify a range like this: RANGE: [from, to, nr]");
            long nr(roundl(range[2]));
            if (nr > MAX_VARIATION_NUMBER) THROW(invalid_input, "You are trying too many Variations, please reconsider.");
            double from(range[0]), to(range[1]);
            values.clear();
            if (to < from) std::swap(to, from);
            values.push_back(from);
            if (to == from) {variation_map->insert(std::make_pair(name, values)); return;}
            double step = (to - from)/(nr-2);
            for (int i = 1; i < nr-1; ++i){
                values.push_back(from + i*step);
            }
            values.push_back(to);
        }
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
    void Variations::FindDependentVertices() {
        msg_Debugging() << "Finding Depending Coulings..." << std::endl;
        // for all the varied parameters
        for (std::string var_name : variables){
            // empty set init
            dependent_vertices.insert(std::make_pair(var_name, new std::set<MODEL::Single_Vertex*>()));
            // go through vertices TODO is this passed by reference???
            std::vector<MODEL::Single_Vertex>* p_vertices = MODEL::s_model->Vertices_Pointer();
            for (std::vector<MODEL::Single_Vertex>::iterator it_v = p_vertices->begin(); it_v != p_vertices->end(); ++it_v) {
                MODEL::Single_Vertex* v = it_v.base(); 
                if (v->DependsOn(var_name)) dependent_vertices[var_name]->emplace(v);
            }
        }
    }

    void Variations::StoreNominal() {
        msg_Debugging() << "Store Nominal Paramter values... " << std::endl;
        std::vector<double_t> values = {};
        for (std::string name : variables){
            if (MODEL::s_model->Constants()->count(name) != 1) return;
            values.push_back(MODEL::s_model->Constants()->at(name));
        }
        nominal = VariationKey(variables, values);
    }

    // singleton
    Variations* Variations::variations = nullptr;
}

