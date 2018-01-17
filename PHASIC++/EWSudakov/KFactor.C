#include "PHASIC++/EWSudakov/KFactor.H"
#include "PHASIC++/EWSudakov/Sudakov.H"
#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;


Sudakov_KFactor::Sudakov_KFactor(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args),
  m_runas(MODEL::One_Running_AlphaS(MODEL::One_Running_AlphaS(p_proc->Integrator()->ISR()->PDF(0)))),
  m_runaqed(MODEL::Running_AlphaQED(1./137.03599976)),p_ews(new Sudakov(p_proc)),
  m_mom(p_proc->Integrator()->Momenta())
{ }
Sudakov_KFactor::~Sudakov_KFactor()
{
  if(p_ews) delete p_ews;
}
double Sudakov_KFactor::KFactor(const int mode)
{
  m_mom.clear();
  double res(0.), pref(1.);
  m_mom = p_proc->Integrator()->Momenta();
  res = pref*p_ews->EWSudakov(m_mom);
  return m_weight=res;
}

double Sudakov_KFactor::KFactor(const ATOOLS::NLO_subevt &evt)
{
  return m_weight=1.;
}

DECLARE_GETTER(Sudakov_KFactor,"EWSudakov",
               KFactor_Setter_Base, KFactor_Setter_Arguments);

KFactor_Setter_Base *ATOOLS::Getter
<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new Sudakov_KFactor(args);
}

void ATOOLS::Getter<KFactor_Setter_Base,KFactor_Setter_Arguments,Sudakov_KFactor>::
PrintInfo(std::ostream &str, const size_t width) const
{
  str << "EW Sudakov K-Factor is implemented in ref ... \n";
}
