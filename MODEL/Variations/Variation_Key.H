#ifndef MODEL_Variations__Variation_Key_H
#define MODEL_Variations__Variation_Key_H

#include "ATOOLS/Org/Exception.H"
#include <math.h>

namespace MODEL {
    namespace VARIATIONS {
        /*
        Type for specifying a single parameter variation using parameter names and corresponding values
        */
        class VariationKey {
            public:
                /* Constructors */
                VariationKey() : names({}), values({}) {id = "empty key";}                          // default empty constructor @param None @returns VariationKey
                VariationKey(std::string var_name, double_t var_value)                              // constructor from one parameter @param var_name string, name of the parameter @param var_value double, value of the parameter @returns VariationKey
                    : names({var_name}), values({var_value}) {UpdateIdentifier();}                  
                VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values); // constructor from multiple parameters @param var_names vector<string>, names of the parameters @param var_values vector<double>, values of the parameters @warning vectors must have same length @returns VariationKey 
                
                /* Getters */
                inline size_t Size() const {return names.size();}                                   // get the size of this key @param None @returns size_t
                inline double Value(size_t idx) const {return values[idx];}                         // get the parameter value at idx @param idx size_t @returns double
                inline const std::vector<double_t>& Values() const {return values;}                 // get the list of parameter values @param None @returns vector<double>
                inline const std::string& Name(size_t idx) const {return names[idx];}               // get parameter name at idx @param idx size_t @returns string
                inline const std::vector<std::string>& Names() const {return names;}                // get the list of parameter names @param None @returns vector<string>
                inline const std::string& Identifier() const {return id;}                           // gets identifier @param None @returns string

                /* Add Parameters */
                void Add(std::string var_name, double_t value);                                     // add a single name and value to this key @param var_name string @param value double @return None
                void Add(std::vector<std::string> var_names, std::vector<double_t> var_values);     // add a list of names and values to this key @param var_names vector<string> @param var_values vector<double> @return None
                
                /* Compare */
                bool const operator< (const VariationKey& other) const;                             // compare with other key by comparing names and then values @param other VariationKey @return bool

            protected:
                /* Attributes */
                std::vector<std::string> names;     // parameter names, one for each value
                std::vector<double> values;         // parameter values, one for each name
                std::string id;                     // identifier
                
                /* Internal Functions */
                void UpdateIdentifier();            // update the identifier string @param None @return None
        };

        /* 
        Output to stream
        @param ss ostream
        @param key VariationKey
        @returns ostream
        */
        inline std::ostream& operator<<(std::ostream& ss, VariationKey key) {ss << key.Identifier(); return ss;}

        /* 
        Add two keys by concatenating their lists of parameters and values.  
        @param k1 VariationKey
        @param k2 VariationKey
        @returns VariationKey
        */
        VariationKey operator+(const VariationKey& k1, const VariationKey& k2);     
    }
}

#endif