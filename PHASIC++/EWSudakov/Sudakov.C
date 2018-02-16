#include "PHASIC++/EWSudakov/Sudakov.H"

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/Coefficient_Checker.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <cassert>

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Sudakov::Sudakov(Process_Base* proc):
  p_proc{ proc },
  m_ampls{ p_proc },
  m_comixinterface{ p_proc, m_ampls },
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) },
  m_check{ Default_Reader().Get<bool>("CHECK_EWSUDAKOV", false) },
  m_coeffs{ {"L", {}}, {"lZ", {}} }
{
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);
  m_SU2rotatedspinampls.clear();
  m_spinampls.clear();
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.Unrotated());
  CalculateSpinAmplitudeCoeffs();
  if (m_check) {
    Coefficient_Checker checker(p_proc->Name());
    if (!checker.CheckCoeffs(m_coeffs, m_spinampls[0])) {
      THROW(fatal_error, "EWSudakov coeffs for this process are not equal to"
                         " the results in hep-ph/0010201.");
    }
  }
  // TODO: dress with logs, calculate factor for squared amplitude
  return 1.0;
}

void Sudakov::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto spinamplnum = ampls.size();
  for (auto& kv : m_coeffs) {
    kv.second.clear();
    kv.second.resize(spinamplnum);
  }
  for (size_t i{ 0 }; i < spinamplnum; ++i) {
    const auto value = ampls.Get(i);
    const auto spincombination = ampls.GetSpinCombination(i);
    if (spincombination.size() != m_ampls.NumberOfLegs())
      THROW(fatal_error, "Inconsistent state");
    if (value == 0.0)
      continue;
    m_coeffs["L"][i] = LsCoeff(value, spincombination, i);
    m_coeffs["lZ"][i] = lsZCoeff(value, spincombination, i);
    // TODO: add other coefficients (remember to init m_coeffs in ctor)
  }
  for (const auto& coeff : m_coeffs["L"])
    DEBUG_VAR(coeff);
}

Complex Sudakov::LsCoeff(Complex amplvalue,
                         std::vector<int> spincombination,
                         size_t spinidx)
{
  Complex coeff{ 0.0 };
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.Unrotated().Leg(i)->Flav() };
    coeff -= DiagonalCew(flav, spincombination[i]) / 2.0;
    if (flav.IsVector() && flav.Charge() == 0 && spincombination[i] != 2) {
      // special case of neutral gauge bosons, they mix and hence non-diagonal
      // terms appear, cf. e.g. eq. (6.30)
      const auto from = flav.Kfcode();
      const auto to = (from == kf_photon) ? kf_Z : kf_photon;
      const auto prefactor = -NondiagonalCew() / 2.0;
      const auto amplratio = 0.0;
      auto amplit = m_SU2rotatedspinampls.find(i);
      auto& legpermutation = m_ampls.LegPermutation(i);
      if (amplit == m_SU2rotatedspinampls.end()) {
        auto& rotatedampl(m_ampls.Rotated(i));
        m_comixinterface.FillSpinAmplitudes(m_SU2rotatedspinampls[i],
                                            rotatedampl);
        amplit = m_SU2rotatedspinampls.find(i);
      }
      std::vector<int> rotatedspincombination;
      for (const auto& idx : legpermutation)
        rotatedspincombination.push_back(spincombination[idx]);
      const auto rotated = amplit->second[0].Get(rotatedspincombination);
      const auto unrotated = amplvalue;
      assert(unrotated != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
      // TODO: minus sign not understood here!
      coeff -= prefactor * rotated / unrotated;
    }
  }
  return coeff;
}

Complex Sudakov::lsZCoeff(Complex amplvalue,
                          std::vector<int> spincombination,
                          size_t spinidx)
{
  Complex coeff{ 0.0 };
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.Unrotated().Leg(i)->Flav() };
    coeff += IZ2(flav, spincombination[i]) * std::log(1.0/m_cw2);
  }
  return coeff;
}

double Sudakov::DiagonalCew(const Flavour& flav, int pol) const
{
  // pol is either chirality or polarisation:
  // 0: + (right-handed or transverse polarisation)
  // 1: - (left-handed or transverse polarisation)
  // 2: 0 (longitudinal polarisation)
  // NOTE: for longitudinal bosons, use the Goldstone equivalence theorem
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
    if (pol == 2)
      return CewLefthandedLepton;
    else
      return 2/m_sw2;
  } else if (flav.IsBoson() && flav.Charge() == 0) {
    if (pol == 2) {
      assert(!flav.IsPhoton());
      return CewLefthandedLepton;
    } else if (flav.IsPhoton()) {
      return 2.0;
    } else {
      return 2.0 * m_cw2/m_sw2;
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double Sudakov::NondiagonalCew() const
{
  return -2.0 * m_cw/m_sw;
}

double Sudakov::IZ2(const Flavour& flav, int pol) const
{
  static auto IZ2LefthandedLepton = std::pow(m_cw2 - m_sw2, 2) / (4*m_sw2*m_cw2);
  static auto IZ2Neutrino = 1 / (4*m_sw2*m_cw2);
  if (flav.IsLepton()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        THROW(fatal_error, "Right-handed neutrino are not supported");
      return m_sw2/m_cw2;
    } else {
      if (flav.IsUptype())
        return IZ2Neutrino;
      else
        return IZ2LefthandedLepton;
    }
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        return 4*m_sw2 / (9*m_cw2);
      else
        return 1*m_sw2 / (9*m_cw2);
    } else {
      if (flav.IsUptype())
        return std::pow(3*m_cw2 - m_sw2, 2) / (36*m_sw2*m_cw2);
      else
        return std::pow(3*m_cw2 + m_sw2, 2) / (36*m_sw2*m_cw2);
    }
  } else if (flav.Kfcode() == kf_Wplus) {
    if (pol == 2)
      return IZ2LefthandedLepton;
    else
      return m_cw2/m_sw2;
  } else if (flav.IsBoson() && flav.Charge() == 0) {
    if (pol == 2) {
      assert(!flav.IsPhoton());
      return IZ2Neutrino;
    } else {
      return 0.0;
    }
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}
