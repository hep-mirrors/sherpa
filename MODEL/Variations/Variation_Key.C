#include "MODEL/Variations/Variation_Key.H"

namespace MODELVARIATIONS {
    /*
    Constructor from 1 parameter with value
    */
    VariationKey::VariationKey(std::string var_name, double_t var_value) {
        names = {var_name}; values = {var_value};
        UpdateIdentifier();
    }

    /*
    Constructor from multiple parameters with one value each
    */
    VariationKey::VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values){
        if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
        names = var_names; values = var_values;
        UpdateIdentifier();
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
    Add another parameter with value to the existing variation and return the new Key
    */
    VariationKey VariationKey::AddVariation(std::string var_name, double_t var_value) {
        std::vector<std::string> new_names = std::vector<std::string>(names);
        std::vector<double_t> new_values = std::vector<double_t>(values);
        new_names.push_back(var_name); new_values.push_back(var_value);
        return VariationKey(new_names, new_values);
    }
}