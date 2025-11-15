#ifndef AddOns_UFOVariations_Combinations_H
#define AddOns_UFOVariations_Combinations_H

#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <algorithm>
#include "AddOns/UFOVariations/UFO_Variation_Key.H"

namespace UFOVariations {
    /*
    Class for combining multiple parameter variations together from a read in map
    */
    class Combinations {
        public:
            Combinations(std::map<std::string, std::vector<double_t>>  *variations_map);
            ~Combinations();
            void AddAllCombinations(std::vector<VariationKey> *variations);

        private:
            std::map<std::string, std::vector<double_t>> *p_variations_map;
            void combine_names(std::vector<std::vector<std::string>> *name_combinations);
            void combine_values(std::vector<std::string> combination, std::vector<std::vector<double_t>> *result_lists);
    };
};

#endif