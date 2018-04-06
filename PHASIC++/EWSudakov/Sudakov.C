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
  m_activecoeffs{ EWSudakov_Log_Type::Ls, EWSudakov_Log_Type::lZ, EWSudakov_Log_Type::lSSC },
  m_ampls{ p_proc, m_activecoeffs },
  m_comixinterface{ p_proc, m_ampls },
  m_sw2{ MODEL::s_model->ComplexConstant("csin2_thetaW").real() },
  m_cw2{ 1.0 - m_sw2 },
  m_sw{ sqrt(m_sw2) },
  m_cw{ sqrt(m_cw2) },
  m_mw2{ sqr(s_kftable[kf_Wplus]->m_mass) },
  m_mz2{ sqr(s_kftable[kf_Z]->m_mass) },
  m_check{ Default_Reader().Get<bool>("CHECK_EWSUDAKOV", false) }
{
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);
  m_lsczspinampls.clear();
  m_sscwspinampls.clear();
  m_spinampls.clear();
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.Unrotated());
  CalculateSpinAmplitudeCoeffs();

  /*
    // pref = alpha/4 pi is added in KFactor method.

    Born_EW = |M0 + alpha * delta * M0 / 4 pi |^2 
            = |M0|^2 |1 + alpha * delta / 4 pi|^2 
            = |M0|^2 (1 + 2 * alpha * Re(delta) / 4 pi) + O(alpha^2)
    This function gives 2 Re(delta), the rest is handled by KFactor.C
  */

  const auto born { p_proc->Get<COMIX::Single_Process>()->Getdxs_beforeKFactor() };
  return 2.*deltaEW((mom[0]+mom[1]).Abs2()).real()/born;
}

Complex Sudakov::deltaEW(const double s)
{
  // TODO: include m_angularcoeffs
  /*
    Define the various logs -> store them in 
    dict with the same names as the coefficients
    then sum all contributions * that log.
   */
  const auto L  = (s>m_mw2)?sqr(std::log(s/m_mw2)):0.;
  const auto lZ = (s>m_mw2)?std::log(s/m_mz2):0.;

  // for now this is a bit of an overkill, but useful
  // if we add more logs...
  std::map<Coeff_Map_Key, double> logs;
  logs[{EWSudakov_Log_Type::Ls, {}}] = L;
  logs[{EWSudakov_Log_Type::lZ, {}}] = lZ;

  Complex res{ 0. };

  for (const auto& coeffkv : m_coeffs) {
    for (const auto c_val : coeffkv.second) {
      res += c_val.first * logs[coeffkv.first];
    }
  }

  return res;
}

void Sudakov::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto spinamplnum = ampls.size();
  m_coeffs.clear();
  for (const auto& key : m_activecoeffs) {
    switch (key) {
      case EWSudakov_Log_Type::Ls:
      case EWSudakov_Log_Type::lZ:
        m_coeffs[{key, {}}].resize(spinamplnum);
        break;
      case EWSudakov_Log_Type::lSSC:
        for (size_t k{ 0 }; k < ampls.GetSpinCombination(0).size(); ++k)
          for (size_t l{ 0 }; l < k; ++l)
            m_coeffs[{key, {k, l}}].resize(spinamplnum);
        break;
    }
  }
  for (size_t i{ 0 }; i < spinamplnum; ++i) {
    const auto value = ampls.Get(i);
    const auto spincombination = ampls.GetSpinCombination(i);
    if (spincombination.size() != m_ampls.NumberOfLegs())
      THROW(fatal_error, "Inconsistent state");
    if (value == 0.0)
      continue;
    for (const auto& key : m_activecoeffs) {
      switch (key) {
        case EWSudakov_Log_Type::Ls:
          m_coeffs[{key, {}}][i] = LsCoeff(value, spincombination, i);
          break;
        case EWSudakov_Log_Type::lZ:
          m_coeffs[{key, {}}][i] = lsZCoeff(value, spincombination, i);
          break;
        case EWSudakov_Log_Type::lSSC:
          for (size_t k{ 0 }; k < spincombination.size(); ++k) {
            for (size_t l{ 0 }; l < k; ++l) {
              const auto angularkey
                = Coeff_Map_Key{EWSudakov_Log_Type::lSSC, {k, l}};
              m_coeffs[angularkey][i]
                = lsLogROverSCoeffs(value, spincombination, i, {k, l});
            }
          }
          break;
      }
    }
  }
  if (m_check) {
    Coefficient_Checker checker(p_proc->Name(), m_activecoeffs);
    if (!checker.CheckCoeffs(m_coeffs, m_spinampls[0])) {
      THROW(fatal_error, "EWSudakov coeffs for this process are not equal to"
                         " the results in hep-ph/0010201.");
    }
  }
  for (size_t i{ 0 }; i < spinamplnum; ++i) {
    const auto value = ampls.Get(i);
    if (value == 0.0)
      continue;
    Complex B0i{ value * std::conj(value) };
    for (auto& coeffkv : m_coeffs)
      coeffkv.second[i].first *= B0i;
  }
}

Coeff_Value Sudakov::LsCoeff(Complex amplvalue,
                             std::vector<int> spincombination,
                             size_t spinidx)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.Unrotated().Leg(i)->Flav() };
    const auto diagonal = DiagonalCew(flav, spincombination[i]) / 2.0;
    coeff.first -= diagonal;
    coeff.second -= diagonal;
    if (flav.IsVector() && flav.Charge() == 0 && spincombination[i] != 2) {
      // special case of neutral gauge bosons, they mix and hence non-diagonal
      // terms appear, cf. e.g. eq. (6.30)
      const auto prefactor = -NondiagonalCew() / 2.0;
      auto amplit = m_lsczspinampls.find(i);
      if (amplit == m_lsczspinampls.end()) {
        auto& rotatedampl =
            m_ampls.Rotated(EWSudakov_Amplitude_Type::LSCZ, {i});
        m_comixinterface.FillSpinAmplitudes(m_lsczspinampls[i], rotatedampl);
        amplit = m_lsczspinampls.find(i);
      }
      auto& legpermutation = m_ampls.LegPermutation(
          EWSudakov_Amplitude_Type::LSCZ, {i});
      std::vector<int> rotatedspincombination;
      for (const auto& idx : legpermutation)
        rotatedspincombination.push_back(spincombination[idx]);
      const auto rotated = amplit->second[0].Get(rotatedspincombination);
      const auto unrotated = amplvalue;
      assert(unrotated != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
      coeff.first += prefactor * rotated / unrotated;
      coeff.second -= prefactor * rotated / unrotated;
    }
  }
  return coeff;
}

Coeff_Value Sudakov::lsZCoeff(Complex amplvalue,
                              std::vector<int> spincombination,
                              size_t spinidx)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.Unrotated().Leg(i)->Flav() };
    // 1/m_cw2 = (mZ/mW)^2 !!! 
    const auto contrib = IZ2(flav, spincombination[i]) * std::log(1.0/m_cw2);
    coeff.first += contrib;
    coeff.second += contrib;
  }
  return coeff;
}

Coeff_Value Sudakov::lsLogROverSCoeffs(Complex amplvalue,
                                       std::vector<int> spincombination,
                                       size_t spinidx,
                                       const Two_Leg_Indizes& indizes)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });

  const auto k = indizes[0];
  const auto l = indizes[1];
  auto kflav = m_ampls.Unrotated().Leg(k)->Flav();
  auto lflav = m_ampls.Unrotated().Leg(l)->Flav();

  // add contribution for each vector boson connecting the leg pairs

  // photon
  const auto IAk = -kflav.Charge();
  const auto IAl = -lflav.Charge();
  coeff.first += 2*IAk*IAl;
  coeff.second += 2*IAk*IAl;

  // Z
  const auto IZk = IZ(kflav, spincombination[k]);
  const auto IZl = IZ(lflav, spincombination[l]);
  coeff.first += 2*IZk*IZl;
  coeff.second += 2*IZk*IZl;

  // W
  for (int i{ 0 }; i < 2; ++i) {
    const auto kplus = (i == 0);
    const auto Ik = Ipm(kflav, spincombination[k], kplus);
    const auto Il = Ipm(lflav, spincombination[l], !kplus);
    if (Ik != 0.0 && Il != 0.0) {
      msg_Debugging() << "calc {" << k << ", " << l << "} i=" << i << std::endl;
      // TODO: remove duplication when calculating ampl ratios
      auto amplit = m_sscwspinampls.find(indizes);
      if (amplit == m_sscwspinampls.end()) {
        auto& rotatedampl =
          m_ampls.Rotated(EWSudakov_Amplitude_Type::SSCW, {k, l});
        m_comixinterface.FillSpinAmplitudes(m_sscwspinampls[indizes], rotatedampl);
        amplit = m_sscwspinampls.find(indizes);
      }
      auto& legpermutation = m_ampls.LegPermutation(
          EWSudakov_Amplitude_Type::SSCW, {k, l});
      std::vector<int> rotatedspincombination;
      for (const auto& idx : legpermutation)
        rotatedspincombination.push_back(spincombination[idx]);
      const auto rotated = amplit->second[0].Get(rotatedspincombination);
      const auto unrotated = amplvalue;
      assert(unrotated != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
      // TODO: understand why we need to use abs here when calculating coeffs
      // for ee->mumu, but when we want to calculate coeffs for ee->uu/dd, we
      // can use the (expected) unmodified ratio; is this connected to the
      // extraneous minus sign in Sudakov::LsCoeff?
      //const auto amplratio = std::abs(rotated/unrotated);
      const auto amplratio = rotated/unrotated;
      DEBUG_VAR(amplratio);
      coeff.first += 2*Ik*Il*amplratio;
      const auto amplratio2
        = Complex{ std::abs(amplratio.real()), amplratio.imag() };
      coeff.second += 2*Ik*Il*amplratio2;
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

double Sudakov::IZ(const Flavour& flav, int pol) const
{
  const auto sign = (flav.IsAnti() ? -1 : 1);
  if (flav.IsScalar())
    THROW(not_implemented,
          "non-diagonal Z coupling terms for scalars not implemented");
  if (flav.IsLepton()) {
    if (pol == 0)
      return sign * m_sw/m_cw;
    else
      return sign * (m_sw2 - m_cw2)/(2*m_cw*m_sw);
  } else if (flav.IsQuark()) {  // cf. eq. (B.16)
    if (pol == 0) {
      if (flav.IsUptype())
        return -sign * 2/3.0 * m_sw/m_cw;
      else
        return sign * 1/3.0 * m_sw/m_cw;
    } else {
      if (flav.IsUptype())
        return sign * (3*m_cw2 - m_sw2) / (6*m_sw*m_cw);
      else
        return -sign * (3*m_cw2 + m_sw2) / (6*m_sw*m_cw);
    }
  } else if (flav.Kfcode() == kf_Wplus) {
    if (pol == 2) {
      // TODO: understand extra minus sign here to get same results as in ref
      // with the Coeff Checker
      return -sign * (m_sw2 - m_cw2)/(2*m_cw*m_sw);
    } else {
      return sign * m_cw/m_sw;
    }
  } else {
    MyStrStream s;
    s << "Missing implementation for flavour: " << flav;
    THROW(not_implemented, s.str());
  }
}

double Sudakov::Ipm(const Flavour& flav, int pol, bool isplus) const
{
  if (pol == 0)
    return 0.0;
  if (flav.IsFermion()) {
    const auto isfermionplus = flav.IsUptype();
    if (flav.IsAnti() && (isplus == isfermionplus))
      return -1 / (sqrt(2)*m_sw);
    else if (!flav.IsAnti() && (isplus != isfermionplus))
      return  1 / (sqrt(2)*m_sw);
    return 0.0;
  } else if (flav.Kfcode() == kf_Wplus) {
    // TODO: This is just a placeholder
    return 0.0;
  } else {
    MyStrStream s;
    s << "Missing implementation for flavour: " << flav;
    THROW(not_implemented, s.str());
  }
}

namespace PHASIC {

  std::ostream& operator<<(std::ostream& os, const Coeff_Map_Key& k)
  {
    os << k.first;
    if (!k.second.empty()) {
      os << " { ";
      for (const auto& i : k.second)
        os << i << " ";
      os << "}";
    }
    return os;
  }

}
