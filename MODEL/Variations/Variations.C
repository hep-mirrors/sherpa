#include "MODEL/Variations/Variations.H"

#define MAX_VARIATION_NUMBER 1000

namespace MODELVARIATIONS {
    /*
    Constructs the Variations Object
    Iterates over the "UFO_VARIATIONS" Setting and adds the Variations to its own list
    */
    Variations::Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner){
        if (!vertices_pointer && !constants_poiner) THROW(fatal_error, "no vertices or constants given.")
        p_vertices = vertices_pointer;
        p_constants = constants_poiner;
        // Deal with UFO Param Variations
        // get settings and read them to map
        msg_Out() << std::endl << "Reading the parameter variations..." << std::endl;
        ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
        std::vector<ATOOLS::Scoped_Settings> items = s["UFO_VARIATIONS"].GetItems();
        m_variations = std::map<std::string, std::vector<double_t>>();
        for (ATOOLS::Scoped_Settings item : items) ReadSingleParamVariation(item);
        // map is done here, add variation key elements from combinations
        AddAllCombinations();        
        // check if there are too many variations
        if (Size() >= MAX_VARIATION_NUMBER) THROW(invalid_input, "You are trying too many Variations, please reconsider.");
        if (msg_LevelIsDebugging()) {
            msg_Out() << "All read Variations are: " << std::endl;
            for (VariationKey key : v_variations){
                msg_Out() << "\t" << key << std::endl;
            }
        }
        else {
            msg_Out() << "\t--> Found " << Size() << " variations" << std::endl;
        }
        // Register the dependent Model vertices, check for unknown/unused variations
        FindDependentVertices();
        for (std::string var_name : variables) {
            if (dependent_vertices[var_name]->size() == 0) THROW(normal_exit, var_name + " does not seem to have anything depending on it. Remove Variation or fix this.") 
        }
        // Store the nominal parameter values 
        StoreNominal();
        msg_Out() << "UFO Variations read." << std::endl << std::endl;
    }

    /*
    Destructor, does something
    */
    Variations::~Variations(){
        for (auto set : dependent_vertices) {
            delete set.second;
        }
    }

    /*
    read a single param variation (one variable, mulitple values) and put it in a map
    */
    void Variations::ReadSingleParamVariation(ATOOLS::Scoped_Settings& s){
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
            if (to == from) {m_variations.insert(std::make_pair(name, values)); return;}
            double step = (to - from)/(nr-2);
            for (int i = 1; i < nr-1; ++i){
                values.push_back(from + i*step);
            }
            values.push_back(to);
        }
        m_variations.insert(std::make_pair(name, values));
    }

    /*
    find and register the dependencies of the couplings of all Model vertices
    */
    void Variations::FindDependentVertices() {
        msg_Debugging() << "Finding Depending Coulings..." << std::endl;
        // for all the varied parameters
        for (std::string var_name : variables){
            // empty set init
            dependent_vertices.insert(std::make_pair(var_name, new std::set<MODEL::Single_Vertex*>()));
            // go through vertices;
            for (std::vector<MODEL::Single_Vertex>::iterator it_v = p_vertices->begin(); it_v != p_vertices->end(); ++it_v) {
                MODEL::Single_Vertex* v = it_v.base(); 
                if (v->DependsOn(var_name)) dependent_vertices[var_name]->emplace(v);
            }
        }
    }

    /*
    Store the nominal parameter values for later reset
    */
    void Variations::StoreNominal() {
        msg_Debugging() << "Store Nominal Paramter values... " << std::endl;
        std::vector<double_t> values = {};
        for (std::string name : variables){
            if (p_constants->count(name) != 1) return;
            values.push_back(p_constants->at(name));
        }
        nominal = VariationKey(variables, values);
    }

    /*
    Adds all possible combinations of parameters and values from the map to the VariationKey vector found at the given pointer.
    */
    void Variations::AddAllCombinations(){
        msg_Debugging() << "Calculate and add all combinations..." << std::endl;
        // save all possible combinations of names to a vector
        std::vector<std::vector<std::string>> name_combinations = {};
        CombineNames(&name_combinations);
        // for each combination of variable names find all the value combinations
        for (std::vector<std::string> combination : name_combinations) {
            // combine them to a list of all possibilities
            std::vector<std::vector<double_t>> combined_val_lists = {};
            CombineValues(combination, &combined_val_lists);
            // now save a key to the passed variations list for each value combination
            for (std::vector<double_t> values : combined_val_lists){
                v_variations.push_back(VariationKey(combination, values));
            }
        }
        msg_Debugging() << "...done" << std::endl;
    }

    /*
    construct all the variable names combinations from a sorted list of names
    */
    void Variations::CombineNames(std::vector<std::vector<std::string>> *name_combinations){
        // get sorted list of all variable names from the map
        std::vector<std::string> sorted_names = {};
        for (auto item : m_variations){
            sorted_names.push_back(item.first);
        }
        std::sort(sorted_names.begin(), sorted_names.end());
        *name_combinations = {};
        for (std::string name : sorted_names)
            name_combinations->push_back({name});
        int i = 0;
        while (i < sorted_names.size()) {
            std::vector<std::vector<std::string>> next_result = std::vector<std::vector<std::string>>(*name_combinations);
            for (std::vector<std::string> list : *name_combinations) {
                if (list.size() == i) {
                    for (std::string name : sorted_names){
                        if (list.back() < name) {
                            std::vector<std::string> new_list = std::vector<std::string>(list);
                            new_list.push_back(name);
                            next_result.push_back(new_list);
                        }
                    }
                }
            }
            *name_combinations = next_result;
            i++;
        }
    }

    /*
    construct all the value combinations from specific lists
    */
    void Variations::CombineValues(std::vector<std::string> combination, std::vector<std::vector<double_t>> *result_lists){
        // find the possible values for each variable name from the map
        std::vector<std::vector<double_t>> lists = {};
        for (std::string var_name : combination) {
            for (auto item : m_variations){
                if (item.first == var_name) lists.push_back(item.second);
            }
        }
        *result_lists = {};
        result_lists->push_back({});
        int i = 0;
        while (i < lists.size()){
            std::vector<std::vector<double_t>> next_result = {};
            std::vector<double_t> list = lists.at(i);
            for (std::vector<double_t> r_list : *result_lists) {
                for (double_t val : list) {
                    std::vector<double_t> new_list = std::vector<double_t>(r_list);
                    new_list.push_back(val);
                    next_result.push_back(new_list);
                }
            }
            *result_lists = next_result;
            i++;
        }
    }
}

