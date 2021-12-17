#include "AddOns/EWSud/KFactor.H"
#include "AddOns/EWSud/Comix_Interface.H"

#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace EWSud;


Sudakov_KFactor::Sudakov_KFactor(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args),
  m_calc{ p_proc }
{
  auto& s = Settings::GetMainSettings();
  m_maxweight = s["EWSUD"]["MAX_KFACTOR"].SetDefault(10.0).Get<double>();
  if(Settings::GetMainSettings()["EWSUDAKOV_MAX_KFACTOR"].IsCustomised()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD:MAX_KFACTOR");
  }
}

double Sudakov_KFactor::KFactor(const int mode)
{
  Calculate();
  Validate();
  return m_weight;
}

double Sudakov_KFactor::KFactor(const ATOOLS::NLO_subevt &evt)
{
  return m_weight = 1.0;
}

void Sudakov_KFactor::CalculateAndFillWeightsMap(Weights_Map& w)
{
  Calculate();
  Validate();
  w["EWSud"]["KFactor"] = m_weight;
  w["EWSud"]["KFactorExp"] = m_expweight;
  for (const auto t : ActiveLogTypes()) {
    w["EWSud"][ToString<EWSudakov_Log_Type>(t)] = 1.0 + m_corrections_map[t];
  }
}

void Sudakov_KFactor::ResetWeightsMap(Weights_Map& w)
{
  w["EWSud"]["KFactor"] = 1.0;
  w["EWSud"]["KFactorExp"] = 1.0;
  for (const auto t : ActiveLogTypes()) {
    w["EWSud"][ToString<EWSudakov_Log_Type>(t)] = 1.0;
  }
}

void Sudakov_KFactor::Calculate()
{
  m_corrections_map = m_calc.CorrectionsMap(p_proc->Integrator()->Momenta());
  m_weight = m_corrections_map.KFactor();
  m_expweight = exp(m_weight - 1.0);
}

void Sudakov_KFactor::Validate()
{
  if (std::abs(m_weight) > m_maxweight) {
    m_weight = 1.0;
  }
  if (std::abs(m_expweight) > m_maxweight) {
    m_expweight = 1.0;
  }
}

DECLARE_GETTER(Sudakov_KFactor,"EWSud",
               KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new Sudakov_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
PrintInfo(std::ostream &str, const size_t width) const
{
  str << "EWSud is implemented in arXiv:2006.14635 and arXiv:2111.13453.\n";
}
