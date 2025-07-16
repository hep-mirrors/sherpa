#ifndef ATOOLS_Phys_Hard_Process_Variation_Generator_Base_H
#define ATOOLS_Phys_Hard_Process_Variation_Generator_Base_H

#include "ATOOLS/Org/Getter_Function.H"

namespace PHASIC { class Process_Base; }

namespace ATOOLS {

  struct Weights_Map;

  struct Hard_Process_Variation_Generator_Arguments {
    PHASIC::Process_Base* p_proc;
  };

  namespace {
    using Args = Hard_Process_Variation_Generator_Arguments;
  }

  class Hard_Process_Variation_Generator_Base {

  public:

    virtual ~Hard_Process_Variation_Generator_Base() = default;

    virtual void GenerateAndFillWeightsMap(Weights_Map&) = 0;
    virtual void ResetWeightsMap(Weights_Map&) = 0;

    static void ShowSyntax(const size_t i);

  };

  using Hard_Process_Variation_Generator_Getter_Function =
    ATOOLS::Getter_Function<Hard_Process_Variation_Generator_Base, Args>;

}

#endif
