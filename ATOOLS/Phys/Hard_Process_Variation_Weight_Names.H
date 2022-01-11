#ifndef ATOOLS_Phys_Hard_Process_Variation_Weight_Names_Base_H
#define ATOOLS_Phys_Hard_Process_Variation_Weight_Names_Base_H

#include "ATOOLS/Org/Getter_Function.H"

namespace ATOOLS {

  struct Hard_Process_Variation_Weight_Names_Arguments {
  };

  namespace {
    using Args = Hard_Process_Variation_Weight_Names_Arguments;
  }

  class Hard_Process_Variation_Weight_Names_Base {

  public:

    virtual ~Hard_Process_Variation_Weight_Names_Base() = default;

    virtual std::string WeightMapKey() const = 0;
    virtual std::string WeightNameForVariation(const std::string&) const = 0;

    static void ShowSyntax(const size_t i);

  };

  using Hard_Process_Variation_Weight_Names_Getter_Function =
    ATOOLS::Getter_Function<Hard_Process_Variation_Weight_Names_Base, Args>;

}

#endif
