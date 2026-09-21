#include "MODEL/Variations/Variation_Key.H"

namespace MODEL {
    namespace VARIATIONS {
        VariationKey::VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values) {
            // check inputs
            if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
            m_names = var_names; m_values = var_values;
            UpdateIdentifier();
        }

        void VariationKey::Add(std::vector<std::string> var_names, std::vector<double_t> var_values){
            if (var_names.size() != var_values.size()) THROW(invalid_input, "input lists are incompatible");
            m_names.insert(m_names.end(), var_names.begin(), var_names.end());
            m_values.insert(m_values.end(), var_values.begin(), var_values.end());
            UpdateIdentifier();
            return;
        }
        
        void VariationKey::Add(std::string var_name, double_t var_value){
            m_names.push_back(var_name);
            m_values.push_back(var_value);
            UpdateIdentifier();
            return;
        }

        void VariationKey::UpdateIdentifier() {
            std::stringstream ss;
            for (int i = 0; i < Size(); i++) {ss << m_names[i] << "-" << m_values[i] << "-.";}
            m_id = ss.str();
        }

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

        VariationKey operator+(const VariationKey& k1, const VariationKey& k2){
            VariationKey result (k1);
            result.Add(k2.Names(), k2.Values());
            return result;
        }
    }
}