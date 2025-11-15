#include "AddOns/UFOVariations/Combinations.H"

namespace UFOVariations {
    /*
    Constructor of the Combinations Object, needs a map of the variable names to their possible values (as a pointer)
    */
    Combinations::Combinations(std::map<std::string, std::vector<double_t>> *variations_map){
        msg_Debugging() << "Built combinations Object" << std::endl;
        p_variations_map = variations_map;
    }

    /*
    Destructor, does nothing
    */
    Combinations::~Combinations() {}
    
    /*
    Adds all possible combinations of parameters and values from the map to the VariationKey vector found at the given pointer.
    */
    void Combinations::AddAllCombinations(std::vector<VariationKey> *variations){
        msg_Debugging() << "Calculate and add all combinations..." << std::endl;
        // save all possible combinations of names to a vector
        std::vector<std::vector<std::string>> name_combinations = {};
        combine_names(&name_combinations);
        // for each combination of variable names find all the value combinations
        for (std::vector<std::string> combination : name_combinations) {
            // combine them to a list of all possibilities
            std::vector<std::vector<double_t>> combined_val_lists = {};
            combine_values(combination, &combined_val_lists);
            // now save a key to the passed variations list for each value combination
            for (std::vector<double_t> values : combined_val_lists){
                variations->push_back(VariationKey(combination, values));
            }
        }
        msg_Debugging() << "...done" << std::endl;
    }

    /*
    construct all the variable names combinations from a sorted list of names
    */
    void Combinations::combine_names(std::vector<std::vector<std::string>> *name_combinations){
        // get sorted list of all variable names from the map
        std::vector<std::string> sorted_names = {};
        for (auto item : *p_variations_map){
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
    void Combinations::combine_values(std::vector<std::string> combination, std::vector<std::vector<double_t>> *result_lists){
        // find the possible values for each variable name from the map
        std::vector<std::vector<double_t>> lists = {};
        for (std::string var_name : combination) {
            for (auto item : *p_variations_map){
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
};