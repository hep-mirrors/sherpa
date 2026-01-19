#ifndef MODEL_Variations_Variations_H
#define MODEL_Variations_Variations_H

#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/Variations/Variation_Key.H"
#include <vector>
#include <map>
#include <set>

namespace MODELVARIATIONS{
    class Variations {
        public:            
            Variations(std::vector<MODEL::Single_Vertex>* vertices_pointer, std::map<std::string, double_t>* constants_poiner);
            ~Variations();
    
            inline size_t Size() const {return v_variations.size();}
            inline VariationKey Nominal() const {return nominal;}
            inline bool IsOkay() const {return okay;}
            inline std::vector<VariationKey> GetVariations() const {return v_variations;};
            inline const std::vector<std::string>& VariableNames() const {return variables;}
            inline std::set<MODEL::Single_Vertex*>* Dependents(std::string var_name) const {if (dependent_vertices.count(var_name) != 1) return new std::set<MODEL::Single_Vertex*>(); return dependent_vertices.at(var_name);};

        protected:
            std::vector<MODEL::Single_Vertex>* p_vertices = nullptr;
            std::map<std::string,double> * p_constants = nullptr;
            ATOOLS::Settings& s = ATOOLS::Settings::GetMainSettings();
            VariationKey nominal = VariationKey::DefaultKey();
            std::vector<std::string> variables;
            std::map<std::string, std::vector<double_t>> m_variations;
            std::vector<VariationKey> v_variations;
            std::map<std::string, std::set<MODEL::Single_Vertex*>*> dependent_vertices;
            bool okay = true;

            void FindDependentVertices();
            void CheckForUnusedVertices();

            void StoreNominal();

            void ReadVariations();
            void ReadSingleParamVariation(std::string parameter);

            void CombineNames(std::vector<std::vector<std::string>> *name_combinations);
            void CombineValues(std::vector<std::string> combination, std::vector<std::vector<double_t>> *result_lists);
            void AddAllCombinations();
            void CheckVariationNumber();
    };
};

#endif