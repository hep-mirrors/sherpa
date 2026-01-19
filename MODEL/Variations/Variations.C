#include "MODEL/Variations/Variations.H"

#define MAX_VARIATION_NUMBER 1000

namespace MODELVARIATIONS {
    /*
    Constructs the Variations Object
    Iterates over the "UFO_VARIATIONS" Setting and adds the Variations to its own list
    */
    Variations::Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner){
        // check arguments and init attributes
        if (!vertices_pointer && !constants_poiner) THROW(fatal_error, "no vertices or constants given by the model")
        p_vertices = vertices_pointer;
        p_constants = constants_poiner;       
        m_variations = std::map<std::string, std::vector<double_t>>();
        variables = std::vector<std::string>();
        // read settings
        msg_Out() << std::endl << "\x1b[34mReading the model parameter variations...\x1b[0m" << std::endl;
        ReadVariations();
        // Register the dependent Model vertices, check for unknown/unused variations
        FindDependentVertices();
        CheckForUnusedVertices();
        // add variation key elements from combinations
        AddAllCombinations();         
        // Store the nominal parameter values 
        StoreNominal();
        // check if there are too many variations
        CheckVariationNumber();
        msg_Out() << "\x1b[34mModel Variations read.\x1b[0m" << std::endl << std::endl;
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
    void Variations::ReadSingleParamVariation(std::string parameter){
        ATOOLS::Scoped_Settings ss = s["MODEL_VARIATIONS"][parameter].SetDefault(-1.);
        std::vector<double_t> values = std::vector<double_t>();
        if (ss.IsMap()){
            DEBUG_INFO("was map");
            if (!(ss["From"].SetDefault(-1.).IsScalar() && ss["To"].SetDefault(-1.).IsScalar() && ss["Step"].SetDefault(-1.).IsScalar())) {
                msg_Out() << "\x1b[31m\tVariations of " << parameter << ": range not formatted properly. Ignoring it...\x1b[0m" << std::endl;
                variables.pop_back();
                return;
            }
            double from = ss["From"].GetScalar<double_t>();
            double to = ss["To"].GetScalar<double_t>();
            double step = ss["Step"].GetScalar<double_t>();
            if (!(from > 0 && to > 0 && step > 0 && from < to)) {
                msg_Out() << "\x1b[31m\tRange for " << parameter << " is not specified correctly. Ignoring it...\x1b[0m" << std::endl;
                variables.pop_back();
                return;
            }
            double val = from;
            // this is supposed to fixed FP errors
            while (val - to <= step/2) {
                values.push_back(val);
                val += step;
            }
        }
        else if (ss.IsList()) {
            DEBUG_INFO("was list");
            values = ss.GetVector<double_t>();
        }
        else if (ss.IsScalar()) {
            DEBUG_INFO("was scalar");
            values.push_back(ss.GetScalar<double_t>());
        }
        if (values.empty() || values.back() < 0) {
            msg_Out() << "\x1b[31m\tVariations of " << parameter << " not formatted properly. Ignoring it...\x1b[0m" << std::endl;
            variables.pop_back();
            return;
        }
        m_variations.insert(std::make_pair(parameter, values));
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
    go over the parameters and remove those who dont have dependents
    */
    void Variations::CheckForUnusedVertices() {
        for (auto var_it=variables.begin(); var_it != variables.end(); var_it++) {
            if (dependent_vertices[*var_it]->size() == 0){
                msg_Out() << "\t\x1b[31m" << *var_it << " does not seem to have anything depending on it. Ignoring it...\x1b[0m" << std::endl;
                variables.erase(var_it);
                var_it--;
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
    Read the parameters from settings
    */
    void Variations::ReadVariations() {
        for (const std::string& parameter : s["MODEL_VARIATIONS"].GetKeys()) {
            DEBUG_INFO("Reading variations of " + parameter);
            variables.push_back(parameter);
            ReadSingleParamVariation(parameter);
        }
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
        std::vector<std::string> sorted_names (variables);
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

    /*
    Checks the Size of the Variations to make sure everything is alright
    */
    void Variations::CheckVariationNumber(){
        if (Size() >= MAX_VARIATION_NUMBER){
            msg_Out() << "\x1b[31m\tYou are trying too many Variations, please reconsider. Ignoring variations...\x1b[0m" << std::endl;
            okay = false;
        }
        else if (Size() == 0) {
            msg_Out() << "\x1b[31m\tNo useful variations are found, ignoring them...\x1b[0m" << std::endl;
            okay = false;
        }
        else {
            // Output
            if (msg_LevelIsDebugging()) {
                msg_Out() << "All read Variations are: " << std::endl;
                for (VariationKey key : v_variations){
                    msg_Out() << "\t" << key << std::endl;
                }
            }
            else {
                msg_Out() << "\x1b[34m\t--> Found " << Size() << " variations.\x1b[0m" << std::endl;
            }
        }
    }
}

