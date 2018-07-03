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
  m_activecoeffs{
    EWSudakov_Log_Type::Ls,
    EWSudakov_Log_Type::lZ,
    EWSudakov_Log_Type::lSSC,
    EWSudakov_Log_Type::lC,
    EWSudakov_Log_Type::lYuk
  },
  m_ampls{ p_proc, m_activecoeffs },
  m_comixinterface{ p_proc, m_ampls },
  m_mw2{ sqr(s_kftable[kf_Wplus]->m_mass) },
  m_mz2{ sqr(s_kftable[kf_Z]->m_mass) },
  m_check{ Default_Reader().Get<bool>("CHECK_EWSUDAKOV", false) }
{
}

double Sudakov::EWSudakov(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);
  if (!IsInHighEnergyLimit())
    return 1.0;

  m_lsczspinampls.clear();
  m_sscwspinampls.clear();
  m_spinampls.clear();
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.BaseAmplitude());
  CalculateSpinAmplitudeCoeffs();
  THROW(normal_exit, "Finish.");

  /*
    // pref = alpha/4 pi is added in KFactor method.

    Born_EW = |M0 + alpha * delta * M0 / 4 pi |^2 
            = |M0|^2 |1 + alpha * delta / 4 pi|^2 
            = |M0|^2 (1 + 2 * alpha * Re(delta) / 4 pi) + O(alpha^2)
    This function gives 2 Re(delta), the rest is handled by KFactor.C
  */

  const auto born = p_proc->Get<COMIX::Single_Process>()->Getdxs_beforeKFactor();
  return 2.*deltaEW((mom[0]+mom[1]).Abs2()).real()/born;
}

bool Sudakov::IsInHighEnergyLimit()
{
  DEBUG_FUNC("");
  static const auto threshold = 1e2;
  const auto s = std::abs(m_ampls.MandelstamS());
  const auto t = std::abs(m_ampls.MandelstamT());
  const auto u = std::abs(m_ampls.MandelstamU());
  DEBUG_VAR(t/s);
  DEBUG_VAR(u/s);

  const auto LSC = sqr(std::log(s/m_mw2));
  const auto SSCt = std::abs(2 * std::log(s/m_mw2) * std::log(std::abs(t)/s));
  const auto SSCu = std::abs(2 * std::log(s/m_mw2) * std::log(std::abs(u)/s));

  msg_Debugging() << "s = " << s << ", t = " << t << ", u = " << u << "\n";
  msg_Debugging() << "log2(s/mW) = " << LSC << "\n";
  msg_Debugging() << "2*log(s/mW)*log(t/s) = " << SSCt << "\n";
  msg_Debugging() << "2*log(s/mW)*log(u/s) = " << SSCu << "\n";

  if (LSC < threshold) {
    msg_Debugging() << "event LSC too small\n";
    return false;
  }
  if ((SSCt < threshold) && (SSCu < threshold)) {
    msg_Debugging() << "event SSC too small\n";
    return false;
  }

  return true;
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
      case EWSudakov_Log_Type::lC:
      case EWSudakov_Log_Type::lYuk:
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
          m_coeffs[{key, {}}][i] = LsCoeff(value, spincombination);
          break;
        case EWSudakov_Log_Type::lZ:
          m_coeffs[{key, {}}][i] = lsZCoeff(value, spincombination);
          break;
        case EWSudakov_Log_Type::lSSC:
          if ((m_ampls.BaseAmplitude().Leg(2)->Flav().Kfcode() == kf_Z
               && spincombination[2] == 2)
              ||
              (m_ampls.BaseAmplitude().Leg(3)->Flav().Kfcode() == kf_Z
               && spincombination[3] == 2)) {
            msg_Error() << "EWSudakov WARNING: omitting SSC coeff calc for ";
            msg_Error() << "an eeZZ (with Z longitudinal) , for now due to ";
            msg_Error() << "missing implementations\n";
            break;
          }
          for (size_t k{ 0 }; k < spincombination.size(); ++k) {
            for (size_t l{ 0 }; l < k; ++l) {
              // s-channel-related loops will have vanishing log coeffs
              if (k == 1 && l == 0)
                continue;
              if (spincombination.size() == 4 && k == 3 && l == 2)
                continue;
              const auto angularkey
                = Coeff_Map_Key{EWSudakov_Log_Type::lSSC, {k, l}};
              m_coeffs[angularkey][i]
                = lsLogROverSCoeffs(value, spincombination, {k, l});
            }
          }
          break;
        case EWSudakov_Log_Type::lC:
          m_coeffs[{key, {}}][i] = lsCCoeff(value, spincombination);
          break;
        case EWSudakov_Log_Type::lYuk:
          m_coeffs[{key, {}}][i] = lsYukCoeff(value, spincombination);
          break;
      }
    }
  }
  if (m_check) {
    Coefficient_Checker checker(p_proc->Name(), m_activecoeffs);
    Mandelstam_Variables mandelstam {
      m_ampls.MandelstamS(),
      m_ampls.MandelstamT(),
      m_ampls.MandelstamU() };
    if (!checker.CheckCoeffs(m_coeffs, m_spinampls[0], mandelstam)) {
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
                             std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    const auto diagonal
      = -m_ewgroupconsts.DiagonalCew(flav, spincombination[i]) / 2.0;
    coeff.first += diagonal;
    coeff.second += diagonal;
    if (flav.IsVector() && flav.Charge() == 0 && spincombination[i] != 2) {
      const kf_code newkf = (flav.Kfcode() == kf_Z) ? kf_photon : kf_Z;
      // special case of neutral gauge bosons, they mix and hence non-diagonal
      // terms appear, cf. e.g. eq. (6.30)
      const auto prefactor = -m_ewgroupconsts.NondiagonalCew() / 2.0;
      auto amplit = m_lsczspinampls.find(i);
      if (amplit == m_lsczspinampls.end()) {
        auto& transformedampl
          = m_ampls.SU2TransformedAmplitude({std::make_pair(i, newkf)});
        m_comixinterface.FillSpinAmplitudes(m_lsczspinampls[i],
                                            transformedampl);
        amplit = m_lsczspinampls.find(i);
      }
      auto& legpermutation = m_ampls.LegPermutation({std::make_pair(i, newkf)});
      std::vector<int> transformedspincombination;
      for (const auto& idx : legpermutation)
        transformedspincombination.push_back(spincombination[idx]);
      auto transformed = amplit->second[0].Get(transformedspincombination);
      const auto base = amplvalue;
      assert(base != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
      auto amplratio = transformed/base;
      coeff.first += prefactor * amplratio;
      coeff.second -= prefactor * amplratio;
    }
  }
  return coeff;
}

Coeff_Value Sudakov::lsZCoeff(Complex amplvalue,
                              std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    // 1/m_cw2 = (mZ/mW)^2 !!! 
    const auto contrib
      = m_ewgroupconsts.IZ2(flav, spincombination[i])
        * std::log(1.0/m_ewgroupconsts.m_cw2);
    coeff.first += contrib;
    coeff.second += contrib;
  }
  return coeff;
}

Coeff_Value Sudakov::lsLogROverSCoeffs(Complex amplvalue,
                                       std::vector<int> spincombination,
                                       const Two_Leg_Indizes& indizes)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });

  const auto k = indizes[0];
  const auto l = indizes[1];
  auto kflav = m_ampls.BaseAmplitude().Leg(k)->Flav();
  auto lflav = m_ampls.BaseAmplitude().Leg(l)->Flav();

  // add contribution for each vector boson connecting the leg pairs

  // photon
  const auto IAk = -kflav.Charge();
  const auto IAl = -lflav.Charge();
  coeff.first += 2*IAk*IAl;
  coeff.second += 2*IAk*IAl;

  // Z
  const auto IZk = m_ewgroupconsts.IZ(kflav, spincombination[k]);
  const auto IZl = m_ewgroupconsts.IZ(lflav, spincombination[l]);
  coeff.first += 2*IZk*IZl;
  coeff.second += 2*IZk*IZl;

  // W
  for (int i{ 0 }; i < 2; ++i) {
    const auto kplus = (i == 0);
    const auto kcouplings
      = m_ewgroupconsts.Ipm(kflav, spincombination[k], kplus);
    const auto lcouplings
      = m_ewgroupconsts.Ipm(lflav, spincombination[l], !kplus);

    for (const auto kcoupling : kcouplings) {
      for (const auto lcoupling : lcouplings) {
        // TODO: remove code duplication when calculating ampl ratios
        const Leg_Set key{ {k, kcoupling.first}, {l, lcoupling.first} };
        auto amplit = m_sscwspinampls.find(key);
        if (amplit == m_sscwspinampls.end()) {
          auto& transformedampl = m_ampls.SU2TransformedAmplitude(key);
          m_comixinterface.FillSpinAmplitudes(m_sscwspinampls[key],
                                              transformedampl);
          amplit = m_sscwspinampls.find(key);
        }

        // correct spin index when a longitudinal vector boson is replaced with
        // a scalar using the Goldstone boson equivalence theorem
        std::vector<int> goldstonespincombination;
        for (size_t i{ 0 }; i < spincombination.size(); ++i) {
          auto lambda = spincombination[i];
          if (lambda == 2) {
            if (i == k && kflav.IsVector() && kcoupling.first != kf_Z) {
              lambda = 0;
            } else if (i == l && lflav.IsVector() && lcoupling.first != kf_Z) {
              lambda = 0;
            }
          }
          goldstonespincombination.push_back(lambda);
        }

        auto& legpermutation = m_ampls.LegPermutation(key);
        std::vector<int> transformedspincombination;
        for (const auto& idx : legpermutation) {
          transformedspincombination.push_back(goldstonespincombination[idx]);
        }
        auto transformed = amplit->second[0].Get(transformedspincombination);
        const auto base = amplvalue;
        assert(base != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
        auto amplratio = transformed/base;
        auto contribution = 2.0*kcoupling.second*lcoupling.second*amplratio;
        const auto amplratio2 = -amplratio;
        auto contribution2 = 2.0*kcoupling.second*lcoupling.second*amplratio2;

        coeff.first += contribution;
        coeff.second += contribution2;
      }
    }
  }
  return coeff;
}

Coeff_Value Sudakov::lsCCoeff(Complex amplvalue,
                              std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i {0}; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    auto contrib = 0.0;
    if (flav.IsFermion()) {
      contrib = 3.0/2.0 * m_ewgroupconsts.DiagonalCew(flav, spincombination[i]);
    } else if (flav.Kfcode() == kf_Wplus && spincombination[i] != 2) {
      contrib = m_ewgroupconsts.Bew(flav, spincombination[i]) / 2.0;
    }
    coeff.first += contrib;
    coeff.second += contrib;
  }
  return coeff;
}


Coeff_Value Sudakov::lsYukCoeff(Complex amplvalue,
                                std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i {0}; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    if (flav.Kfcode() != kf_t && flav.Kfcode() != kf_b)
      continue;
    auto contrib = sqr(flav.Mass()/Flavour{kf_Wplus}.Mass());
    if (spincombination[i] == 0)
      contrib *= 2.0;
    else
      contrib
        += sqr(flav.IsoWeakPartner().Mass()/Flavour{kf_Wplus}.Mass());
    contrib *= -1.0/(8.0*m_ewgroupconsts.m_sw2);
    coeff.first += contrib;
    coeff.second += contrib;
  }
  return coeff;
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
