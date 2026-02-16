#ifndef MODEL_Variations_Variations_H
#define MODEL_Variations_Variations_H

#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Variations/Variation_Key.H"
#include <vector>
#include <map>
#include <set>

namespace MODEL{
    namespace VARIATIONS {
        /*
        Class for reading and storing the model parameters variations
        */
        class Variations {
            public:
                /*
                Constructs the Variations Object. Iterates over the "MODEL_VARIATIONS" settings and adds the variations to its own list.
                @param vertices_pointer pointer to a vector of MODEL::Single_Vertex, usually just the vertices owned by the model.
                @param constants_poiner pointer to a map of string to real constants, usually just the constants owned by the model
                @return Variations
                */
                Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner);
                ~Variations();                                                                      // Destructor, does something @param None @return delete
                inline size_t Size() const {return variations_list.size();}                         // size of variations @param None @return size_t
                inline VariationKey Nominal() const {return nominal;}                               // nominal variation @param None @return VariationKey
                inline bool IsOkay() const {return okay;}                                           // check this object @param None @return bool
                inline std::vector<VariationKey> GetVariations() const {return variations_list;}    // get all variations @param None @return vector<VariationKey>
                inline std::set<MODEL::Single_Vertex*>* Dependents(std::string var_name) const      // gets the dependents on a parameter @param var_name name of the parameter @return set<Single_Vertex*>
                    {if (dependent_vertices.count(var_name) != 1) return new std::set<MODEL::Single_Vertex*>(); return dependent_vertices.at(var_name);} 
                        

            protected:
                /* Attributes */
                std::vector<VariationKey> variations_list;                                          // list of all variations
                VariationKey nominal;                                                               // nominal variation
                std::map<std::string, std::set<MODEL::Single_Vertex*>*> dependent_vertices;         // map of dependent vertices for each parameter
                bool okay = true;                                                                   // status of this object

                /* model information */
                std::vector<MODEL::Single_Vertex>* p_vertices = nullptr;                            // all model vertices
                std::map<std::string,double> * p_constants = nullptr;                               // all model constants

                /* Settings */
                ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();                          // main settings

                /* calculation buffers */
                std::vector<std::string> variable_names;                                            // names of all parameters
                std::map<std::string, std::vector<double_t>> variations_map;                        // map from parameter names to their values
                std::vector<std::vector<std::string>> correlated_variables;                         // lists of correlated variables

                /* internal functions for convenience */
                void FindDependentVertices();                           // find and register the dependencies of the couplings of all model vertices @param None @return None
                void CheckParameters();                                 // go over the parameters and remove those who dont have dependents @param None @return None
                void StoreNominal();                                    // store the nominal parameter values for reset @param None @return None
                void ReadVariations();                                  // read the parameters from settings @param None @return None
                void ReadSingleParamVariation(std::string parameter);   // read a single variation (one parameter, mulitple values) and put it in a map @param parameter string with name of the parameter @return None
                void ReadCorrelations();                                // read correlations as lists of parameter names @param None @return None
                void CombineParameters();                               // combine parameters depending on the settings and correlations @param None @return None
                void CheckVariationNumber();                            // checks the size of the variations @param None @return None
        };
    };
};

#endif