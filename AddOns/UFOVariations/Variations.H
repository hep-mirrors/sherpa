#ifndef AddOns_UFOVariations_Variations_H
#define AddOns_UFOVariations_Variations_H

#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Kabbala.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"
#include "AddOns/UFOVariations/UFO_Variation_Key.H"
#include <vector>
#include <map>

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
            inline std::vector<ATOOLS::Kabbala*> Dependents(std::string var_name) {if (dependent_couplings.count(var_name) != 1) return std::vector<ATOOLS::Kabbala*>(); return dependent_couplings[var_name];};

        protected:
            Variations();
            ~Variations();

            static Variations* variations;
            VariationKey nominal = VariationKey::DefaultKey();
            std::vector<std::string> variables;
            std::vector<VariationKey> v_variations;
            std::map<std::string, std::vector<ATOOLS::Kabbala*>> dependent_couplings;

            void FindDependentCouplings();
            void StoreNominal();
            void ReadSingleParamVariation(ATOOLS::Scoped_Settings& s, std::map<std::string, std::vector<double_t>> *variation_map);
    };
};

#endif