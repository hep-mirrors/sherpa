#ifndef ATOOLS_Phys_Hard_Process_Variation_Generator_Base_H
#define ATOOLS_Phys_Hard_Process_Variation_Generator_Base_H

#include "ATOOLS/Org/Getter_Function.H"

namespace PHASIC { class Process_Base; }

namespace ATOOLS {

  struct Weights_Map;

  struct Hard_Process_Variation_Generator_Arguments {
    PHASIC::Process_Base *p_proc;
  };

  class Hard_Process_Variation_Generator_Base {

  public:

    using Args = Hard_Process_Variation_Generator_Arguments;

    using Hard_Process_Variation_Generator_Getter_Function =
      ATOOLS::Getter_Function<Hard_Process_Variation_Generator_Base, Args>;

    Hard_Process_Variation_Generator_Base(const Args&);

    virtual ~Hard_Process_Variation_Generator_Base() = default;

    virtual void GenerateAndFillWeightsMap(Weights_Map&) = 0;
    virtual void ResetWeightsMap(Weights_Map&) = 0;

    std::string WeightNameForVariation(const std::string&);

    static void ShowSyntax(const size_t i);

  protected:

    PHASIC::Process_Base* p_proc;

  };

}

#endif

