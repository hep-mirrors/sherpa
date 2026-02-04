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
        class Variations {
            public:            
                Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner);
                ~Variations();
        
                inline size_t Size() const {return variations_list.size();}
                inline VariationKey Nominal() const {return nominal;}
                inline bool IsOkay() const {return okay;}
                inline std::vector<VariationKey> GetVariations() const {return variations_list;}
                inline std::set<MODEL::Single_Vertex*>* Dependents(std::string var_name) const {if (dependent_vertices.count(var_name) != 1) return new std::set<MODEL::Single_Vertex*>(); return dependent_vertices.at(var_name);}

            protected:
                // model information
                std::vector<MODEL::Single_Vertex>* p_vertices = nullptr;
                std::map<std::string,double> * p_constants = nullptr;
                // settings are needed everywhere
                ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
                // nominal values
                VariationKey nominal;
                // things for calculation
                std::vector<std::string> variable_names;
                std::map<std::string, std::vector<double_t>> variations_map;
                std::map<std::string, std::set<MODEL::Single_Vertex*>*> dependent_vertices;
                std::vector<std::vector<std::string>> correlated_variables;
                bool okay = true;
                // list of variations to be used
                std::vector<VariationKey> variations_list;

                // functions for more readability
                void FindDependentVertices();
                void CheckParameters();
                void StoreNominal();
                void ReadVariations();
                void ReadSingleParamVariation(std::string parameter);
                void ReadCorrelations();
                void CombineParameters();
                void CheckVariationNumber();
        };
    };
};

#endif