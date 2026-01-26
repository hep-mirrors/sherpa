#ifndef MODEL_Variations__Variation_Key_H
#define MODEL_Variations__Variation_Key_H

#include "ATOOLS/Org/Exception.H"
#include <math.h>

namespace MODEL {
    namespace VARIATIONS {
        class VariationKey {
            public:
                // Constructors
                VariationKey(std::string var_name, double_t var_value);
                VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values);
                
                // Getters
                inline size_t Size() const {return names.size();}
                inline double Value(size_t idx) const {return values[idx];}
                inline std::string Name(size_t idx) const {return names[idx];}
                inline const std::vector<std::string>& Names() const {return names;}
                inline std::string Identifier() const {return id;}
                
                // Compare
                bool const operator< (const VariationKey& other) const;

                // static default
                inline static VariationKey DefaultKey() {return VariationKey{"DEFAULT", 0.0};}

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
    }
}

#endif