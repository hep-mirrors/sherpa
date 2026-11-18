#ifndef AddOns_GoSam_Variation_Generator_H
#define AddOns_GoSam_Variation_Generator_H

#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "PHASIC++/Process/Single_Process.H"

#include <map>
#include <string>
#include <vector>

namespace GoSam {

  class Variation_Generator :
    public ATOOLS::Hard_Process_Variation_Generator_Base {

  public:

    Variation_Generator(
        const ATOOLS::Hard_Process_Variation_Generator_Arguments&);

    void GenerateAndFillWeightsMap(ATOOLS::Weights_Map&) override;

    void ResetWeightsMap(ATOOLS::Weights_Map&) override;

  private:

    struct Variation {
      std::string m_name;
      std::map<std::string,double> m_params;
    };

    PHASIC::Single_Process* p_proc;
    std::vector<Variation> m_variations;
    std::map<std::string,double> m_nominal;

  };

}

#endif
