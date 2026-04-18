#include "AddOns/THRew/Variation_Generator.H"

#include <ostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace THRew;

using Base = Hard_Process_Variation_Generator_Base;
using Args = Hard_Process_Variation_Generator_Arguments;

Variation_Generator::Variation_Generator(const Args& args)
    : m_kfactor{KFactor_Setter_Arguments{"THRew", args.p_proc}}
{
}

void Variation_Generator::GenerateAndFillWeightsMap(Weights_Map& wgtmap)
{
  m_kfactor.CalculateAndFillWeightsMap(wgtmap);
}

void Variation_Generator::ResetWeightsMap(Weights_Map& wgtmap)
{
  m_kfactor.ResetWeightsMap(wgtmap);
}

DECLARE_GETTER(Variation_Generator, "THRew", Base, Args);

Base* ATOOLS::Getter<Base, Args, Variation_Generator>::operator()(
    const Args& args) const
{
  return new Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, Variation_Generator>::PrintInfo(
    std::ostream& str, const size_t width) const
{
  str << "THRew: Theory Reweighter.\n";
}
