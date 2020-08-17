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
  m_check = s["CHECK_EWSUDAKOV"].SetDefault(false).Get<bool>();
}

double Sudakov_KFactor::KFactor(const int mode)
{
  const auto level = msg->Level();
  if (m_check) msg->SetLevel(8);
  m_corrections_map = m_calc.CorrectionsMap(p_proc->Integrator()->Momenta());
  m_weight = m_corrections_map.KFactor();
  if (m_check) msg->SetLevel(level);
  if (std::abs(m_weight) > 500) {
    msg_Info() << "WARNING: K factor is really large: " << m_weight << '\n';
  }
  // for now, use the m_ewsudakovkfacdelta of the ME wgt info to propagate the
  // weight information on a coeff type-to-type basis
  m_weight = 1.0;
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
