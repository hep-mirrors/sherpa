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
  m_calc{ p_proc }, m_maxweight(10.), m_expweight(1.), m_write_contribs(false)
{
  auto& s = Settings::GetMainSettings();
  m_maxweight = s["EWSUD"]["MAX_KFACTOR"].SetDefault(10.0).Get<double>();
  if(Settings::GetMainSettings()["EWSUDAKOV_MAX_KFACTOR"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD:MAX_KFACTOR");
  }
  m_write_contribs = s["EWSUD"]["WRITE_CONTRIBS"].SetDefault(false).Get<bool>();
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
  WriteNominal(w);
  WriteContribs(w);
  WriteThresholdVariations(w);
}

void Sudakov_KFactor::WriteNominal(Weights_Map& w)
{
  if (m_calc.NThresholds()>1) return;
  w["EWSud"]["EWNLL"] = m_weight;
  w["EWSud"]["ExpEWNLL"] = m_expweight;
}

void Sudakov_KFactor::WriteContribs(Weights_Map& w)
{
  if (m_write_contribs) {
    for (const auto t : ActiveLogTypes()) {
      w["EWSud"][ToString<EWSudakov_Log_Type>(t)] = 1.0 + m_corrections_map[t];
    }
  }
}

void Sudakov_KFactor::WriteThresholdVariations(Weights_Map& w)
{
  if (m_calc.NThresholds()==1) return;
  DEBUG_FUNC("n_thr = "<<m_calc.NThresholds());
  // fill weights, we know thresholds are ordered ascending
  // still always fill all weight
  double wgt(m_weight), expwgt(m_expweight);
  bool ishel(true);
  for (double thr : m_calc.Thresholds()) {
    if (ishel && !m_calc.IsInHighEnergyLimit(thr)) {
      ishel = false; wgt = 1.; expwgt = 1.;
    }
    msg_Debugging()<<"thr = "<<thr
                   <<", wgt = "<<wgt<<", exp(wgt) = "<<expwgt<<std::endl;
    w["EWSud"]["EWNLL_Thr"+ToString(thr)] = wgt;
    w["EWSud"]["ExpEWNLL_Thr"+ToString(thr)] = expwgt;
  }
}

void Sudakov_KFactor::ResetWeightsMap(Weights_Map& w)
{
  if (m_calc.NThresholds()==1) {
    w["EWSud"]["EWNLL"] = 1.0;
    w["EWSud"]["ExpEWNLL"] = 1.0;
  }
  else {
    for (double thr : m_calc.Thresholds()) {
      w["EWSud"]["EWNLL_Thr"+ToString(thr)] = 1.0;
      w["EWSud"]["ExpEWNLL_Thr"+ToString(thr)] = 1.0;
    }
  }
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
