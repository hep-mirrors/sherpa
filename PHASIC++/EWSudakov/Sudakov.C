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

#include <cassert>

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Sudakov::Sudakov(Process_Base& proc):
  m_proc{ proc },
  m_ampls{ proc },
  m_ci{ m_proc, &m_ampls.Unrotated() },
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) },
  m_check{ Default_Reader().Get<bool>("CHECK_EWSUDAKOV", false) }
{
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);
  m_SU2rotatedspinampls.clear();
  m_spinampls.clear();
  m_ci.FillSpinAmplitudes(m_spinampls, &m_ampls.Unrotated());
  CalculateSpinAmplitudeCoeffs();
  if (m_check && !CheckCoeffs()) {
    THROW(fatal_error, "EWSudakov coeffs for this process are not equal to"
                       " the results in hep-ph/0010201.");
  }
  // TODO: dress with logs, calculate factor for squared amplitude
  return 1.0;
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
  if (spincombination.size() != m_ampls.NumberOfLegs())
    THROW(fatal_error, "Inconsistent state");
  Complex coeff{ 0.0 };
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.Unrotated().Leg(i)->Flav() };
    coeff -= DiagonalCew(flav, spincombination[i]) / 2.0;
    if (flav.IsVector() && flav.Charge() == 0 && spincombination[i] != 2) {
      // special case of neutral gauge bosons, they mix and hence non-diagonal
      // terms appear, cf. e.g. eq. (6.30)
      const auto from = flav.Kfcode();
      const auto to = (from == kf_photon) ? kf_Z : kf_photon;
      const auto prefactor = -NondiagonalCew(to, from) / 2.0;
      const auto amplratio = 0.0;
      auto amplit = m_SU2rotatedspinampls.find(i);
      auto& legpermutation = m_ampls.LegPermutation(i);
      if (amplit == m_SU2rotatedspinampls.end()) {
        // TODO: this should be moved to a higher level, also make sure that
        // the spin amplitudes are re-calculated for each event
        auto& rotatedampl(m_ampls.Rotated(i));
        m_ci.FillSpinAmplitudes(m_SU2rotatedspinampls[i], &rotatedampl);
        amplit = m_SU2rotatedspinampls.find(i);
      }
      std::vector<int> rotatedspincombination;
      for (const auto& idx : legpermutation)
        rotatedspincombination.push_back(spincombination[idx]);
      const auto rotated = amplit->second[0].Get(rotatedspincombination);
      const auto unrotated = ampls.Get(spinidx);
      assert(unrotated != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
      // TODO: minus sign not understood here!
      coeff -= prefactor * rotated / unrotated;
    } 
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
    if (pol == 2)
      // TODO: confirm that this is true even for photons
      return CewLefthandedLepton;
    else if (flav.IsPhoton())
      return 2.0;
    else
      return 2.0 * m_cw2/m_sw2;
  } else {
    THROW(not_implemented, "Missing implementation");
  }
}

double Sudakov::NondiagonalCew(kf_code from, kf_code to) const
{
  assert((from == kf_Z && to == kf_photon) || (from == kf_photon && to == kf_Z));
  return -2.0 * m_cw/m_sw;
}

bool Sudakov::CheckCoeffs()
{
  auto res = true;
  const auto& refs = ReferenceCoeffs();
  for (const auto& helrefpair : refs) {
    const auto& helicities = helrefpair.first;
    const auto idx = m_spinampls[0].GetNumber(helicities);
    const auto coeff = m_coeffs[idx];
    msg_Debugging() << om::red;
    for (const auto& h : helicities)
      msg_Debugging() << h << " ";
    msg_Debugging()
      << "coeff: " << coeff
      << "\t vs \t  reference value: " << helrefpair.second
      << om::reset << std::endl;
    const auto prec = (std::abs(helrefpair.second) < 10.0) ? 1.e-2 : 1.e-1;
    if (IsBad(coeff.real()) || std::abs(coeff.real() - helrefpair.second) > prec) {
      res = false;
    }
  }
  return res;
}

const Sudakov::HelicityCoeffMap& Sudakov::ReferenceCoeffs()
{
  static std::map<std::string, HelicityCoeffMap> coeffs;
  if (coeffs.empty()) {
    auto& mapmm = coeffs["2_2__e-__e+__mu-__mu+"];
    mapmm[{0, 0, 0, 0}] = -2.58;
    mapmm[{1, 1, 0, 0}] = -4.96;
    mapmm[{0, 0, 1, 1}] = -4.96;
    mapmm[{1, 1, 1, 1}] = -7.35;
    auto& mapuu = coeffs["2_2__e-__e+__u__ub"];
    mapuu[{0, 0, 0, 0}] = -1.86;
    mapuu[{1, 1, 0, 0}] = -4.25;
    mapuu[{0, 0, 1, 1}] = -4.68;
    mapuu[{1, 1, 1, 1}] = -7.07;
    auto& mapdd = coeffs["2_2__e-__e+__d__db"];
    mapdd[{0, 0, 0, 0}] = -1.43;
    mapdd[{1, 1, 0, 0}] = -3.82;
    mapdd[{0, 0, 1, 1}] = -4.68;
    mapdd[{1, 1, 1, 1}] = -7.07;
    auto& mapWW = coeffs["2_2__e-__e+__W+__W-"];
    mapWW[{0, 0, 2, 2}] = -4.96;
    mapWW[{1, 1, 2, 2}] = -7.35;
    mapWW[{1, 1, 0, 1}] = -12.6;
    mapWW[{1, 1, 1, 0}] = -12.6;
    auto& mapPP = coeffs["2_2__e-__e+__P__P"];
    mapPP[{0, 0, 0, 1}] = -1.29;
    mapPP[{0, 0, 1, 0}] = -1.29;
    mapPP[{1, 1, 0, 1}] = -8.15;
    mapPP[{1, 1, 1, 0}] = -8.15;
    auto& mapZP = coeffs["2_2__e-__e+__Z__P"];
    mapZP[{0, 0, 0, 1}] = -1.29;
    mapZP[{0, 0, 1, 0}] = -1.29;
    mapZP[{1, 1, 0, 1}] = -12.2;
    mapZP[{1, 1, 1, 0}] = -12.2;
    auto& mapZZ = coeffs["2_2__e-__e+__Z__Z"];
    mapZZ[{0, 0, 0, 1}] = -1.29;
    mapZZ[{0, 0, 1, 0}] = -1.29;
    mapZZ[{1, 1, 0, 1}] = -16.2;
    mapZZ[{1, 1, 1, 0}] = -16.2;
  }

  // check proc name is inside the few we have
  const auto pname(m_proc.Name());
  size_t check_name(0);
  for(auto it = coeffs.begin(); it != coeffs.end(); ++it) {
    if (it->first == pname) {
      check_name = 1;
      break;
    }
  }
  if (!check_name)
    THROW(not_implemented, "No test for proc: " + pname);

  return coeffs[pname];
}
