#include "MODEL/Variations/Variation_Key.H"

namespace MODEL {
    namespace VARIATIONS {
        /*
        Constructor from multiple parameters with one value each
        */
        VariationKey::VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values) {
            // check inputs
            if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
            names = var_names; values = var_values;
            UpdateIdentifier();
        }

        /*
        Add things
        */
        void VariationKey::Add(std::vector<std::string> var_names, std::vector<double_t> var_values){
            if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
            names.insert(names.end(), var_names.begin(), var_names.end());
            values.insert(values.end(), var_values.begin(), var_values.end());
            UpdateIdentifier();
            return;
        }
        
        /*
        Add one thing
        */
        void VariationKey::Add(std::string var_name, double_t var_value){
            names.push_back(var_name);
            values.push_back(var_value);
            UpdateIdentifier();
            return;
        }

        /*
        set id depending on the list of values
        */
        void VariationKey::UpdateIdentifier() {
            std::stringstream ss;
            for (int i = 0; i < Size(); i++) {ss << names[i] << "-" << values[i] << "-.";}
            id = ss.str();
        }

        /*
        Comparison for usage as map key
        */
        bool const VariationKey::operator< (const VariationKey& other) const {
            if (this->Size() == 0) return other.Size() > 0;
            for (int i = 0; i < this->Size(); i++) {
                if (i >= other.Size()) return true;
                if (this->Name(i) == other.Name(i)) {
                    if (this->Value(i) == other.Value(i)) continue;
                    else return this->Value(i) < other.Value(i);
                }
                else return this->Name(i) < other.Name(i);
            }
            return true;
        }

        /*
        Add variations to each other
        */
        VariationKey operator+(const VariationKey& k1, const VariationKey& k2){
            VariationKey result (k1);
            result.Add(k2.Names(), k2.Values());
            return result;
        }
    }
}