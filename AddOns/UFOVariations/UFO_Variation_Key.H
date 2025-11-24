#ifndef AddOns_UFOVariations__UFO_Variation_Key_H
#define AddOns_UFOVariations__UFO_Variation_Key_H

#include <string>
#include <sstream>
#include "ATOOLS/Org/Exception.H"
#include <math.h>

namespace UFOVariations {
    struct VariationKey
        {
        VariationKey(std::string var_name, double_t var_value) {
            names = {var_name}; values = {var_value};
            UpdateIdentifier();
        }
        VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values){
            if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
            names = var_names; values = var_values;
            UpdateIdentifier();
        }
        std::vector<std::string> names;
        std::vector<double> values;
        std::string id;
        size_t Size() const {return names.size();}
        double Value(size_t idx) const {return values[idx];}
        std::string Name(size_t idx) const {return names[idx];}
        const std::vector<std::string>& Names() const {return names;}
        void UpdateIdentifier() {
            std::stringstream ss;
            for (int i = 0; i < Size(); i++) {ss << names[i] << "-" << values[i] << "-.";}
            id = ss.str();
        }
        std::string Identifier() const {return id;}
        bool const operator< (const VariationKey& other) const {
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
        VariationKey AddVariation(std::string var_name, double_t var_value) {
            std::vector<std::string> new_names = std::vector<std::string>(names);
            std::vector<double_t> new_values = std::vector<double_t>(values);
            new_names.push_back(var_name); new_values.push_back(var_value);
            return VariationKey(new_names, new_values);
        }
        static VariationKey DefaultKey() {return VariationKey{"DEFAULT", 0.0};}
    };

    inline std::ostream& operator<<(std::ostream& s, VariationKey key)
    {
        s << key.Identifier();
        return s;
    }
}

#endif