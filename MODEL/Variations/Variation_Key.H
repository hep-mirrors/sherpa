#ifndef MODEL_Variations__Variation_Key_H
#define MODEL_Variations__Variation_Key_H

#include "ATOOLS/Org/Exception.H"
#include <math.h>

namespace MODEL {
    namespace VARIATIONS {
        /*
        Type for specifying a single parameter variation using parameter m_names and corresponding m_values
        */
        class VariationKey {
            public:
                /* Constructors */
                VariationKey() : m_names({}), m_values({}) {m_id = "empty key";}                        // default empty constructor @param None @returns VariationKey
                VariationKey(std::string var_name, double_t var_value)                                  // constructor from one parameter @param var_name string, name of the parameter @param var_value double, value of the parameter @returns VariationKey
                    : m_names({var_name}), m_values({var_value}) {UpdateIdentifier();}                  
                VariationKey(std::vector<std::string> var_names, std::vector<double_t> var_values);     // constructor from multiple parameters @param var_m_names vector<string>, m_names of the parameters @param var_m_values vector<double>, m_values of the parameters @warning vectors must have same length @returns VariationKey 
                
                /* Getters */
                inline size_t Size() const {return m_names.size();}                                   // get the size of this key @param None @returns size_t
                inline double Value(size_t idx) const {return m_values[idx];}                         // get the parameter value at idx @param idx size_t @returns double
                inline const std::vector<double_t>& Values() const {return m_values;}                 // get the list of parameter m_values @param None @returns vector<double>
                inline const std::string& Name(size_t idx) const {return m_names[idx];}               // get parameter name at idx @param idx size_t @returns string
                inline const std::vector<std::string>& Names() const {return m_names;}                // get the list of parameter m_names @param None @returns vector<string>
                inline const std::string& Identifier() const {return m_id;}                           // gets identifier @param None @returns string

                /* Setters */
                inline void SetIdentifier(const std::string& identifier) {m_id = identifier;}

                /* Add Parameters */
                void Add(std::string var_name, double_t value);                                     // add a single name and value to this key @param var_name string @param value double @return None
                void Add(std::vector<std::string> var_names, std::vector<double_t> var_values);     // add a list of m_names and m_values to this key @param var_m_names vector<string> @param var_m_values vector<double> @return None
                
                /* Compare */
                bool const operator< (const VariationKey& other) const;                             // compare with other key by comparing m_names and then m_values @param other VariationKey @return bool

            protected:
                /* Attributes */
                std::vector<std::string> m_names;     // parameter m_names, one for each value
                std::vector<double> m_values;         // parameter m_values, one for each name
                std::string m_id;                     // identifier
                
                /* Internal Functions */
                void UpdateIdentifier();              // update the identifier string @param None @return None
        };

        /* 
        Output to stream
        @param ss ostream
        @param key VariationKey
        @returns ostream
        */
        inline std::ostream& operator<<(std::ostream& ss, VariationKey key) {ss << key.Identifier(); return ss;}

        /* 
        Add two keys by concatenating their lists of parameters and m_values.  
        @param k1 VariationKey
        @param k2 VariationKey
        @returns VariationKey
        */
        VariationKey operator+(const VariationKey& k1, const VariationKey& k2);     
    }
}

#endif