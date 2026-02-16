#include "MODEL/Variations/Variations.H"

#define MAX_VARIATION_NUMBER 1000

namespace MODEL {
    namespace VARIATIONS {
        Variations::Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner){
            // check arguments and init attributes
            if (!vertices_pointer && !constants_poiner) THROW(fatal_error, "no vertices or constants given by the model")
            p_vertices = vertices_pointer;
            p_constants = constants_poiner;
            // read settings
            msg_Out() << std::endl << "Reading the model parameter variations..." << std::endl;
            ReadVariations();
            // Register the dependent Model vertices, check for unknown/unused variations
            FindDependentVertices();
            CheckParameters();
            // find wanted correlations, using MODEL_VARIATIONS_CORRELATE
            ReadCorrelations();
            // add variation key elements from combinations
            CombineParameters();
            // Store the nominal parameter values 
            StoreNominal();
            // check if there are too many variations
            CheckVariationNumber();
            msg_Out() << "Model Variations read." << std::endl << std::endl;
        }

        Variations::~Variations(){
            for (auto set : dependent_vertices) {
                delete set.second;
            }
        }
        
        void Variations::ReadCorrelations(){
            if (!s["MODEL_VARIATIONS_CORRELATE"].SetDefault("None").IsList()) {
                if (s["MODEL_VARIATIONS_CORRELATE"].Get<std::string>() != "None")
                    msg_Out() << "\x1b[31m\tParameter correlations of not formatted properly. Ignoring them...\x1b[0m" << std::endl;
                return;
            }
            for (ATOOLS::Scoped_Settings& ss : s["MODEL_VARIATIONS_CORRELATE"].GetItems()){
                if (!ss.IsList()) {
                    msg_Out() << "\x1b[31m\tParameter correlations of not formatted properly. Ignoring them...\x1b[0m" << std::endl;
                    return;
                }
                std::vector<std::string> variables = ss.GetVector<std::string>();
                int length = -1;
                for (const std::string& var : variables){
                    // check if the variable names exist in variations
                    bool exists = false;
                    for (const std::string& var2 : variable_names){
                        if (var == var2) {exists = true; break;};
                    }
                    if (!exists) {
                        msg_Out() << "\x1b[31m\tParameter " << var << " specified for correlations does not exist. Ignoring correlations...\x1b[0m" << std::endl;
                        return;
                    }
                    // check if correlation is specfied twice (illegal)
                    for (const std::vector<std::string>& list : correlated_variables){
                        for (const std::string& var2 : list){
                            if (var == var2) {exists = false; break;};
                        }
                    }
                    if (!exists) {
                        msg_Out() << "\x1b[31m\tParameter " << var << " was specified twice for correlations, which does not make sense. Ignoring correlations...\x1b[0m" << std::endl;
                        return;
                    }
                    // check for same length > 0
                    if (length == -1) length = variations_map[var].size();
                    if (length != variations_map[var].size() || length == 0) {
                        msg_Out() << "\x1b[31m\tParameter " << var << " has a different number of variations than the others (or is zero) for correlations, which does not make sense. Ignoring correlations...\x1b[0m" << std::endl;
                        return;
                    }
                }   
                correlated_variables.push_back(variables);
                msg_Out() << std::endl << "\x1b[34m\tFound correlation: "<< variables <<"\x1b[0m" << std::endl;
            }
        }

        void Variations::ReadSingleParamVariation(std::string parameter){
            ATOOLS::Scoped_Settings ss = s["MODEL_VARIATIONS"][parameter].SetDefault(-1.);
            std::vector<double_t> values = std::vector<double_t>();
            if (ss.IsMap()){
                DEBUG_INFO("was map");
                if (!(ss["Min"].SetDefault(-1.).IsScalar() && ss["Max"].SetDefault(-1.).IsScalar() && ss["Step"].SetDefault(-1.).IsScalar())) {
                    msg_Out() << "\x1b[31m\tVariations of " << parameter << ": range not formatted properly. Ignoring it...\x1b[0m" << std::endl;
                    variable_names.pop_back();
                    return;
                }
                double from = ss["Min"].GetScalar<double_t>();
                double to = ss["Max"].GetScalar<double_t>();
                double step = ss["Step"].GetScalar<double_t>();
                if (!(from > 0 && to > 0 && step > 0 && from < to)) {
                    msg_Out() << "\x1b[31m\tRange for " << parameter << " is not specified correctly. Ignoring it...\x1b[0m" << std::endl;
                    variable_names.pop_back();
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
                variable_names.pop_back();
                return;
            }
            variations_map.insert(std::make_pair(parameter, values));
        }

        void Variations::FindDependentVertices() {
            msg_Debugging() << "Finding Depending Coulings..." << std::endl;
            // for all the varied parameters
            for (std::string var_name : variable_names){
                // empty set init
                dependent_vertices.insert(std::make_pair(var_name, new std::set<MODEL::Single_Vertex*>()));
                // go through vertices;
                for (std::vector<MODEL::Single_Vertex>::iterator it_v = p_vertices->begin(); it_v != p_vertices->end(); ++it_v) {
                    MODEL::Single_Vertex* v = it_v.base(); 
                    if (v->DependsOn(var_name)) dependent_vertices[var_name]->emplace(v);
                }
            }
        }

        void Variations::CheckParameters() {
            for (auto var_it=variable_names.begin(); var_it != variable_names.end(); var_it++) {
                if (dependent_vertices[*var_it]->size() == 0){
                    msg_Out() << "\t\x1b[31m" << *var_it << " does not seem to have anything depending on it. Ignoring it...\x1b[0m" << std::endl;
                    variable_names.erase(var_it);
                    var_it--;
                }
                if (p_constants->find(*var_it) == p_constants->end()){
                    msg_Out() << "\t\x1b[31m" << *var_it << " is not found in the real model constants. Ignoring it...\x1b[0m" << std::endl;
                    variable_names.erase(var_it);
                    var_it--;
                }
            }
        }

        void Variations::StoreNominal() {
            msg_Debugging() << "Store Nominal Paramter values... " << std::endl;
            std::vector<double_t> values = {};
            for (std::string name : variable_names){
                if (p_constants->count(name) != 1) return;
                values.push_back(p_constants->at(name));
            }
            nominal = VariationKey(variable_names, values);
        }

        void Variations::ReadVariations() {
            for (const std::string& parameter : s["MODEL_VARIATIONS"].GetKeys()) {
                DEBUG_INFO("Reading variations of " + parameter);
                variable_names.push_back(parameter);
                ReadSingleParamVariation(parameter);
            }
        }

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
                    for (VariationKey key : variations_list){
                        msg_Out() << "\t" << key << std::endl;
                    }
                }
                else {
                    msg_Out() << "\t--> Found " << Size() << " variations." << std::endl;
                }
            }
        }

        void Variations::CombineParameters(){
            std::vector<std::vector<VariationKey>> individual_variations_list;
            std::vector<VariationKey> individual_variations;
            // deal with the correlations first, reading is dealt with earlier
            for (const std::vector<std::string>& list : correlated_variables){
                individual_variations.clear();
                int length = variations_map[list.front()].size();
                for (int i = 0; i < length; i++){
                    std::vector<double_t> values;
                    for (const std::string& variable : list){
                        values.push_back(variations_map[variable][i]);
                    }
                    individual_variations.push_back(VariationKey(list, values));
                }
                // add nominal values as a key
                VariationKey key = VariationKey();
                for (const std::string& variable : list) key.Add(variable, p_constants->at(variable));
                individual_variations.push_back(key);
                individual_variations_list.push_back(individual_variations);
            }
            // now push all other uncorrelated parameters to the individuals list
            for (const std::string& variable : variable_names){
                individual_variations.clear();
                // check if it is correlated, if so skip it
                bool skip = false;
                for (const std::vector<std::string>& list : correlated_variables) {
                    for (const std::string& var2 : list) {
                        if (variable == var2) {skip = true; break;}
                    }
                    if (skip) break;
                }
                if (skip) continue;
                // now push to the list
                for (double_t value : variations_map[variable]) individual_variations.push_back(VariationKey(variable, value));
                // add nominal
                individual_variations.push_back(VariationKey(variable, p_constants->at(variable)));
                individual_variations_list.push_back(individual_variations);
            }
            // now we are ready to combine the things on the indivduals list
            // first get mode, 0 for no combining, 1 for combining all (default), 2 for one each combination (ignores all correlations and values)
            int mode = 1;
            if (s["MODEL_VARIATIONS_COMBINE"].SetDefault(1.).IsScalar()) mode = s["MODEL_VARIATIONS_COMBINE"].SetDefault(1.).Get<int>();
            // maybe more modes in the future
            switch (mode) {
                case 0: 
                    // no combination, just put all keys into the variations list
                    // dont add the nominals though, which are the last key
                    for (const std::vector<VariationKey>& keys : individual_variations_list)
                        for (std::vector<VariationKey>::const_iterator it = keys.begin(); it != keys.end() - 1; it++)
                            variations_list.push_back(*(it));
                    break;
                case 2:
                    // one each (meaning all zero except 1) for EFT scans
                    // ignore the correlations and the given values
                    // just use nominal value for 1 variable and 0 for the others
                    for (const std::string& var_name1 : variable_names){
                        VariationKey key = VariationKey(var_name1, p_constants->at(var_name1));
                        for (const std::string& var_name2 : variable_names)
                            if (var_name1 != var_name2) key.Add(var_name2, 0.);
                        variations_list.push_back(key);
                    }
                    break;
                default:
                    // combine all the elements from the individual lists
                    // to do this nominal value was added to each list
                    std::vector<VariationKey> new_list;
                    variations_list.push_back(VariationKey());
                    for (int i = 0; i < individual_variations_list.size(); i++){
                        for (const VariationKey& key : individual_variations_list[i]) {
                            for (const VariationKey& old_key : variations_list) {
                                new_list.push_back(key + old_key);
                            }
                        }
                        variations_list = new_list;
                        new_list.clear();
                    }
                    // remove the nominal
                    variations_list.pop_back();
                    break;
            }
        }
    }
}