#include "PHASIC++/EWSudakov/KFactor.H"
#include "PHASIC++/EWSudakov/Comix_Interface.H"

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


Sudakov_KFactor::Sudakov_KFactor(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args),
  m_calc{ p_proc }
{
  auto& s = Settings::GetMainSettings();
}

double Sudakov_KFactor::KFactor(const int mode)
{
  const auto level = msg->Level();
  m_corrections_map = m_calc.CorrectionsMap(p_proc->Integrator()->Momenta());
  m_weight = m_corrections_map.KFactor();
  if (std::abs(m_weight) > 10) {
    msg_Error() << "KFactor from EWSud is large: " << m_weight << " -> ignore\n"
                << m_corrections_map << '\n';
    for (auto& kv : m_corrections_map) {
      kv.second = 0.0;
    }
    m_weight = 1.0;
  }
  return m_weight;
}

double Sudakov_KFactor::KFactor(const ATOOLS::NLO_subevt &evt)
{
  return m_weight = 1.0;
}

DECLARE_GETTER(Sudakov_KFactor,"EWSudakov",
               KFactor_Setter_Base,KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new Sudakov_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
PrintInfo(std::ostream &str, const size_t width) const
{
  str << "EW Sudakov K-Factor is implemented in ref ... \n";
}
