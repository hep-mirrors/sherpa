#ifndef AddOns_UFOVariations_Variations_H
#define AddOns_UFOVariations_Variations_H

#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Kabbala.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "AddOns/UFOVariations/UFO_Variation_Key.H"
#include <vector>
#include <map>
#include <set>

namespace UFOVariations{
    class Variations {
        public:
            static Variations* Get();
            std::vector<VariationKey> GetVariations() const {return v_variations;};
            static size_t Count() {if (variations) return variations->Size(); else return 0;}
            Variations(Variations &other) = delete;
            void operator=(const Variations &) = delete;
            static void Init();
            size_t Size() const {return v_variations.size();}
            VariationKey Nominal() const {return nominal;}
            const std::vector<std::string>& VariableNames() const {return variables;}
            inline std::set<MODEL::Single_Vertex*>* Dependents(std::string var_name) const {if (dependent_vertices.count(var_name) != 1) return new std::set<MODEL::Single_Vertex*>(); return dependent_vertices.at(var_name);};

        protected:
            Variations();
            ~Variations();

            static Variations* variations;
            VariationKey nominal = VariationKey::DefaultKey();
            std::vector<std::string> variables;
            std::vector<VariationKey> v_variations;
            std::map<std::string, std::set<MODEL::Single_Vertex*>*> dependent_vertices;

            void FindDependentVertices();
            void StoreNominal();
            void ReadSingleParamVariation(ATOOLS::Scoped_Settings& s, std::map<std::string, std::vector<double_t>> *variation_map);
    };
};

#endif