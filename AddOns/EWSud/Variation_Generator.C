#include "AddOns/EWSud/Variation_Generator.H"

#include <ostream>

using namespace PHASIC;
using namespace ATOOLS;
using namespace EWSud;

using Base = Hard_Process_Variation_Generator_Base;
using Args = Hard_Process_Variation_Generator_Arguments;

Variation_Generator::Variation_Generator(const Args& args):
  m_kfactor{KFactor_Setter_Arguments{"EWSud", args.p_proc}}
{
  auto s = Settings::GetMainSettings()["EWSUD"];
  m_ewsudakov_rs = s["RS"].SetDefault(true).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_RS"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: RS");
  }
}

void Variation_Generator::GenerateAndFillWeightsMap(Weights_Map& wgtmap)
{
  if (m_kfactor.Process()->GetSubevtList() == nullptr || m_ewsudakov_rs) {
    m_kfactor.CalculateAndFillWeightsMap(wgtmap);
  } else {
    ResetWeightsMap(wgtmap);
  }
}

void Variation_Generator::ResetWeightsMap(Weights_Map& wgtmap)
{
  m_kfactor.ResetWeightsMap(wgtmap);
}

DECLARE_GETTER(Variation_Generator, "EWSud", Base, Args);

Base* ATOOLS::Getter<Base, Args, Variation_Generator>::
operator()(const Args& args) const
{
  return new Variation_Generator(args);
}

void ATOOLS::Getter<Base, Args, Variation_Generator>::
PrintInfo(std::ostream& str, const size_t width) const
{ 
  str << "EWSud is implemented in arXiv:2006.14635 and arXiv:2111.13453.\n";
}

