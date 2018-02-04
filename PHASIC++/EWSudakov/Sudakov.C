#include "PHASIC++/EWSudakov/Sudakov.H"
#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

// TODO: in the final implementation we should read this from the SIN2THETAW
// setting instead, for now this value is useful for comparisons to
// Denner:2000jv
constexpr auto sinw2 = 0.22356;
constexpr auto cosw2 = 1 - sinw2;

Sudakov::Sudakov(Process_Base& proc):
  m_proc{ proc },
  p_ampl{ Sudakov::CreateAmplitude(m_proc) },
  m_ci{ m_proc, *p_ampl }
{ }

/// effective electroweak Casimir operator, cf. eq. (B.10)
double Sudakov::Cew(const Flavour& flav, int pol)
{
  // pol is either chirality or polarisation:
  // 0: + (right-handed or transverse polarisation)
  // 1: - (left-handed or transverse polarisation)
  // 2: 0 (longitudinal polarisation)
  if (flav.IsLepton()) {
    if (pol == 0)
      return 1/cosw2;
    else
      return (1 + 2*cosw2) / (4*sinw2*cosw2);
  } else if (flav.IsQuark()) {
    if (pol == 1) {
      return (sinw2 + 27*cosw2) / (36*sinw2*cosw2);
    } else {
      if (flav.IsUptype())
        return 4 / (9*cosw2);
      else
        return 1 / (9*cosw2);
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

ATOOLS::Cluster_Amplitude* Sudakov::CreateAmplitude(Process_Base& proc)

{
  ATOOLS::Cluster_Amplitude* ampl{ Cluster_Amplitude::New() };
  ampl->SetNIn(proc.NIn());
  ampl->SetOrderQCD(proc.MaxOrder(0));
  for (size_t i(1);i<proc.MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+proc.MaxOrder(i));
  for(int i(0);i<proc.NIn()+proc.NOut();++i)
    if (i<proc.NIn()) ampl->CreateLeg(Vec4D(),proc.Flavours()[i].Bar());
    else ampl->CreateLeg(Vec4D(),proc.Flavours()[i]);
  ampl->SetProc(&proc);
  ampl->SetProcs(proc.AllProcs());
  return ampl;
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  for (size_t i{ 0 }; i < mom.size(); ++i) {
    p_ampl->Leg(i)->SetMom(mom[i]);
  }
  m_ci.FillSpinAmplitudes(m_spinampls);
  const auto DLCoeff = CalculateDoubleLogarithmicCoefficient();
  DEBUG_VAR(DLCoeff);
  return 1.0;
}

double Sudakov::CalculateDoubleLogarithmicCoefficient() const
{
  DEBUG_FUNC("");
  std::vector<double> coeffs(m_spinampls[0].size());
  for (size_t i{ 0 }; i < m_spinampls[0].size(); ++i) {
    const auto value = m_spinampls[0].Get(i);
    if (value == 0.0)
      continue;
    const auto spincombination = m_spinampls[0].GetSpinCombination(i);
    if (spincombination.size() != p_ampl->Legs().size())
      THROW(fatal_error, "Inconsistent state");
    for (size_t j{ 0 }; j < spincombination.size(); ++j) {
      const Flavour flav{ p_ampl->Leg(j)->Flav() };
      if (flav.IsBoson() && flav.Charge() == 0) {
        THROW(not_implemented, "Missing implementation");
      } else {
        coeffs[i] -= Cew(p_ampl->Leg(j)->Flav(), spincombination[j]) / 2.0;
      }
    }
  }
  DEBUG_VAR(coeffs);
  return 0.0;
}
