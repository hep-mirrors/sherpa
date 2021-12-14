#include "AddOns/EWSudakov/Variation_Weight_Names.H"

#include <ostream>

using namespace ATOOLS;
using namespace EWSud;

using Base = Hard_Process_Variation_Weight_Names_Base;
using Args = Hard_Process_Variation_Weight_Names_Arguments;

Variation_Weight_Names::Variation_Weight_Names(const Args& args)
{
}

std::string Variation_Weight_Names::WeightMapKey() const
{
  return "EWSud";
}

std::string
Variation_Weight_Names::WeightNameForVariation(const std::string& varname) const
{
  return WeightMapKey() + "_" + varname;
}

DECLARE_GETTER(Variation_Weight_Names, "EWSud", Base, Args);

Base* ATOOLS::Getter<Base, Args, Variation_Weight_Names>::
operator()(const Args& args) const
{
  return new Variation_Weight_Names(args);
}

void ATOOLS::Getter<Base, Args, Variation_Weight_Names>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
  str << "Weight names for EWSud(akov) variations"; 
}
