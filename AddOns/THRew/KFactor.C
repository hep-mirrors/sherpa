#include "AddOns/THRew/KFactor.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace THRew;

THRew_KFactor::THRew_KFactor(const KFactor_Setter_Arguments& args)
    : KFactor_Setter_Base(args)
{
}

double THRew_KFactor::KFactor(const int mode)
{
  Calculate();
  return m_weight;
}

double THRew_KFactor::KFactor(const ATOOLS::NLO_subevt& evt)
{
  return m_weight = 1.0;
}

void THRew_KFactor::CalculateAndFillWeightsMap(Weights_Map& w)
{
  Calculate();
  w["THRew"]["KFactor"] = m_weight;
}

void THRew_KFactor::ResetWeightsMap(Weights_Map& w)
{
  w["THRew"]["KFactor"] = 1.0;
}

void THRew_KFactor::Calculate()
{
  // TODO: implement theory reweighting calculation
  m_weight = 1.0;
}

DECLARE_GETTER(THRew_KFactor, "THRew", KFactor_Setter_Base,
               KFactor_Setter_Arguments);

KFactor_Setter_Base*
ATOOLS::Getter<KFactor_Setter_Base, KFactor_Setter_Arguments,
               THRew_KFactor>::operator()(const KFactor_Setter_Arguments& args)
    const
{
  return new THRew_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base, KFactor_Setter_Arguments,
                    THRew_KFactor>::PrintInfo(std::ostream& str,
                                              const size_t width) const
{
  str << "THRew: Theory Reweighter K-factor.\n";
}
