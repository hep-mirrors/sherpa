#ifndef MODEL_Variations__Variation_Key_H
#define MODEL_Variations__Variation_Key_H

#include "ATOOLS/Org/Exception.H"
#include <math.h>

namespace MODEL {
    namespace VARIATIONS {
        class VariationKey {
            public:
                // Constructors
                VariationKey() : names({}), values({}) {id = "empty key";}
                VariationKey(std::string var_name, double_t var_value) : names({var_name}), values({var_value}) {UpdateIdentifier();}
                VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values);
                
                // Getters
                inline size_t Size() const {return names.size();}
                inline double Value(size_t idx) const {return values[idx];}
                inline const std::vector<double_t>& Values() const {return values;}
                inline std::string Name(size_t idx) const {return names[idx];}
                inline const std::vector<std::string>& Names() const {return names;}
                inline std::string Identifier() const {return id;}

                // add variations
                void Add(std::string var_name, double_t value);
                void Add(std::vector<std::string> var_names, std::vector<double_t> var_values);
                
                // Compare
                bool const operator< (const VariationKey& other) const;

            protected:
                // attributes
                std::vector<std::string> names;
                std::vector<double> values;
                std::string id;
                
                // Update the id
                void UpdateIdentifier();
        };
        // Output
        inline std::ostream& operator<<(std::ostream& s, VariationKey key) {s << key.Identifier(); return s;}
        // operations
        VariationKey operator+(const VariationKey& k1, const VariationKey& k2);
    }
}

#endif