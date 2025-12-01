#ifndef MODEL_UFO_Variations_H
#define MODEL_UFO_Variations_H

#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "MODEL/UFO/UFO_Variation_Key.H"
#include "MODEL/Main/Model_Base.H"
#include <vector>
#include <map>
#include <set>

namespace UFO{
    class Variations {
        public:            
            Variations();
            ~Variations();
    
            inline size_t Size() const {return v_variations.size();}
            inline VariationKey Nominal() const {return nominal;}
            inline std::vector<VariationKey> GetVariations() const {return v_variations;};
            inline const std::vector<std::string>& VariableNames() const {return variables;}
            inline std::set<MODEL::Single_Vertex*>* Dependents(std::string var_name) const {if (dependent_vertices.count(var_name) != 1) return new std::set<MODEL::Single_Vertex*>(); return dependent_vertices.at(var_name);};

        protected:
            VariationKey nominal = VariationKey::DefaultKey();
            std::vector<std::string> variables;
            std::map<std::string, std::vector<double_t>> m_variations;
            std::vector<VariationKey> v_variations;
            std::map<std::string, std::set<MODEL::Single_Vertex*>*> dependent_vertices;

            void FindDependentVertices();
            void StoreNominal();
            void ReadSingleParamVariation(ATOOLS::Scoped_Settings& s);

            void CombineNames(std::vector<std::vector<std::string>> *name_combinations);
            void CombineValues(std::vector<std::string> combination, std::vector<std::vector<double_t>> *result_lists);
            void AddAllCombinations();
    };
};

#endif