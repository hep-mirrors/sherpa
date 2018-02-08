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

Sudakov::Sudakov(Process_Base& proc):
  m_proc{ proc },
  p_ampl{ Sudakov::CreateAmplitude(m_proc) },
  m_ci{ m_proc, p_ampl },
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) }
{
}

Cluster_Amplitude* Sudakov::CreateAmplitude(Process_Base& proc)
{
  Cluster_Amplitude* ampl{ Cluster_Amplitude::New() };
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

Cluster_Amplitude* Sudakov::CreateSU2RotatedAmplitude(size_t legindex) const
{
  auto* ampl = p_ampl->Copy();
  auto* leg = ampl->Leg(legindex);
  auto flav = leg->Flav();
  Flavour newflav;
  // TODO: generalise to other flavours (preferredly, add an WeakIsoPartner
  // function to the Flavour class
  if (flav.IsPhoton())
    newflav = Flavour{kf_Z};
  else if (flav.Kfcode() == kf_Z)
    newflav = Flavour{kf_photon};
  leg->SetFlav(newflav);
  return ampl;
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  UpdateAmplitude(mom);
  m_SU2rotatedspinampls.clear();  // they will be calculated on demand
  m_ci.FillSpinAmplitudes(m_spinampls, p_ampl);
  CalculateSpinAmplitudeCoeffs();
  // TODO: combine and return coefficients
  return 1.0;
}

void Sudakov::UpdateAmplitude(const ATOOLS::Vec4D_Vector& mom)
{
  for (size_t i{ 0 }; i < mom.size(); ++i) {
    p_ampl->Leg(i)->SetMom(mom[i]);
  }
}

void Sudakov::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto spinamplnum = ampls.size();
  m_coeffs = std::vector<Complex>(spinamplnum, 0.0);
  for (size_t i{ 0 }; i < spinamplnum; ++i) {
    const auto value = ampls.Get(i);
    if (value == 0.0)
      continue;
    m_coeffs[i] = DoubleLogCoeff(m_spinampls[0], i);
    // TODO: add other coefficients
  }
  for (const auto& coeff : m_coeffs)
    DEBUG_VAR(coeff);
}

Complex Sudakov::DoubleLogCoeff(const Spin_Amplitudes& ampls, size_t spinidx)
{
  const auto spincombination = m_spinampls[0].GetSpinCombination(spinidx);
  if (spincombination.size() != p_ampl->Legs().size())
    THROW(fatal_error, "Inconsistent state");
  Complex coeff{ 0.0 };
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ p_ampl->Leg(i)->Flav() };
    if (flav.IsBoson() && flav.Charge() == 0) {
      // mixing between neutral gauge bosons: non-diagonal terms appear
      const auto from = flav.Kfcode();
      coeff -= NondiagonalCew(from, from) / 2.0;
      const auto to = (from == kf_photon) ? kf_Z : kf_photon;
      const auto prefactor = -NondiagonalCew(from, to) / 2.0;
      const auto amplratio = 0.0;
      auto it = m_SU2rotatedspinampls.find(i);
      if (it == m_SU2rotatedspinampls.end()) {
        auto* ampl = CreateSU2RotatedAmplitude(i);
        m_ci.FillSpinAmplitudes(m_SU2rotatedspinampls[i], ampl);
        it = m_SU2rotatedspinampls.find(i);
        ampl->Delete();
      }
      coeff += prefactor * it->second[0].Get(spinidx) / ampls.Get(spinidx);
    } else {
      // only diagonal terms appear
      coeff -= DiagonalCew(flav, spincombination[i]) / 2.0;
    }
  }
  // TODO: multiply with L(s)
  return coeff;
}

double Sudakov::DiagonalCew(const Flavour& flav, int pol) const
{
  // pol is either chirality or polarisation:
  // 0: + (right-handed or transverse polarisation)
  // 1: - (left-handed or transverse polarisation)
  // 2: 0 (longitudinal polarisation)
  static auto CewLefthandedLepton = (1 + 2*m_cw2) / (4*m_sw2*m_cw2);
  if (flav.IsLepton()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return 1/m_cw2;
    } else {
      return CewLefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 1) {
      return (m_sw2 + 27*m_cw2) / (36*m_sw2*m_cw2);
    } else {
      if (flav.IsUptype())
        return 4 / (9*m_cw2);
      else
        return 1 / (9*m_cw2);
    }
  } else if (flav.IsScalar()) {  // cf. eq. (B.18) and (B.16)
    return CewLefthandedLepton;
  } else if (flav.Kfcode() == kf_Wplus) {
    // NOTE: for longitudinal bosons, use the Goldstone equivalence theorem
    if (pol == 2)
      return CewLefthandedLepton;
    else
      return 2/m_sw2;
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double Sudakov::NondiagonalCew(kf_code from, kf_code to) const
{
  if ((from != kf_Z && from != kf_photon) || (to != kf_Z && to != kf_photon))
    THROW(fatal_error, "Only neutral gauge bosons are supported");
  if (from != to)
    return -2.0 * m_cw/m_sw;
  if (from == kf_photon)
    return 2.0;
  if (from == kf_Z)
    return 2.0 * m_cw2/m_sw2;
  THROW(fatal_error, "Logic error");
}
