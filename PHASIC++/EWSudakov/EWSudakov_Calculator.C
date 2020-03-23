#include "PHASIC++/EWSudakov/EWSudakov_Calculator.H"

#include "Coefficient_Checker.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include <cassert>

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Histogram EWSudakov_Calculator::m_kfachisto(0, -5.0, 5.0, 50);
size_t EWSudakov_Calculator::m_numonshellwarning {0};

EWSudakov_Calculator::EWSudakov_Calculator(Process_Base* proc):
  p_proc{ proc },
  m_activecoeffs{
    EWSudakov_Log_Type::Ls,
    EWSudakov_Log_Type::lZ,
    EWSudakov_Log_Type::lSSC,
    EWSudakov_Log_Type::lC,
    EWSudakov_Log_Type::lYuk,
    EWSudakov_Log_Type::lPR
  },
  m_ampls{ p_proc, m_activecoeffs },
  m_comixinterface{ p_proc, m_ampls },
  // TODO: probably we do not want to set up all (SU(2)-transformed) processes
  // again for doing the PR logs, so re-consider what ampls we pass here, such
  // that the HE COMIX interface does not set up processes that we will never
  // need
  m_comixinterface_he{ p_proc, m_ampls }
{
  auto& s = Settings::GetMainSettings();
  m_check = s["CHECK_EWSUDAKOV"].SetDefault(false).Get<bool>();
  m_threshold = s["EWSUDAKOV_THRESHOLD"].SetDefault(5.0).Get<double>();
  m_checkinvariantratios = s["EWSUDAKOV_CHECKINVARIANTRATIOS"].SetDefault(false).Get<bool>();
  s.DeclareVectorSettingsWithEmptyDefault({"EWSUDAKOV_COEFF_REMOVED_LIST"});
  const auto disabled_log_list =
      s["EWSUDAKOV_COEFF_REMOVED_LIST"].GetVector<std::string>();
  for (const auto& l : disabled_log_list) {
    m_activecoeffs.erase(EWSudakovLogTypeFromString(l));
  }
  msg_Out() << "\n ";
  PRINT_INFO("Active EW_Sudakov coefficients : ");
  for (const auto& key : m_activecoeffs)
    msg_Out() << om::red << "\t" << key << om::reset << "\n";
}

EWSudakov_Calculator::~EWSudakov_Calculator()
{
  static bool did_output{false};
  if (!did_output) {
    EWSudakov_Calculator::m_kfachisto.MPISync();
    EWSudakov_Calculator::m_kfachisto.Finalize();
    MyStrStream s;
    s << "kfacs_" << m_threshold;
    EWSudakov_Calculator::m_kfachisto.Output(s.str());
    msg_Error() << "Set " << m_numonshellwarning
                << " amplitudes to 0.0, because there was not enough energy to "
                   "fulfil on-shell conditions\n";
    did_output = true;
  }
}

double EWSudakov_Calculator::KFactor(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  p_proc->GetMEwgtinfo()->m_ewsudakovkfacdelta.clear();
  for (const auto& c : m_activecoeffs)
    p_proc->GetMEwgtinfo()->m_ewsudakovkfacdelta[c] = 0.0;
  m_ampls.UpdateMomenta(mom);
  if (!IsInHighEnergyLimit())
    return 1.0;
  if (p_proc->Integrator()->ColorScheme() == cls::sample) {
    Int_Vector I = p_proc->Integrator()->ColorIntegrator()->I();
    Int_Vector J = p_proc->Integrator()->ColorIntegrator()->J();
    assert(I.size() == J.size());
    assert(I.size() == m_ampls.NumberOfLegs());
    m_ampls.UpdateColors(I, J);
  }
  ClearSpinAmplitudes();
  FillBaseSpinAmplitudes();
  CalculateSpinAmplitudeCoeffs();
  return KFactor();
}

bool EWSudakov_Calculator::IsInHighEnergyLimit()
{
  DEBUG_FUNC("");
  static const auto threshold = sqr(m_threshold) * m_ewgroupconsts.m_mw2;

  const auto s = std::abs(m_ampls.MandelstamS());

  const auto& base_ampl = m_ampls.BaseAmplitude();
  for (size_t i {0}; i < base_ampl.Legs().size(); ++i) {
    for (size_t j {i + 1}; j <  base_ampl.Legs().size(); ++j) {
      const auto sij
        = std::abs((base_ampl.Leg(i)->Mom() + base_ampl.Leg(j)->Mom()).Abs2());
      if(sij < threshold) {
        return false;
      }
      if (m_checkinvariantratios && sij*sij < s * m_ewgroupconsts.m_mw2) {
        return false;
      }
    }
  }
  return true;
}

void EWSudakov_Calculator::ClearSpinAmplitudes()
{
  m_spinampls.clear();
  m_comixinterface.ResetSpinAmplitudeCache();
  m_comixinterface_he.ResetSpinAmplitudeCache();
}

void EWSudakov_Calculator::FillBaseSpinAmplitudes()
{
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.BaseAmplitude());
}

double EWSudakov_Calculator::KFactor()
{
  auto den = m_spinampls[0].SumSquare();
  if (den == 0.0)
    return 1.0;
  const auto s = std::abs(m_ampls.MandelstamS());
  const auto ls = std::log(s/m_ewgroupconsts.m_mw2);

  // pre-calculate the logarithms we need below
  std::map<Coeff_Map_Key, double> logs;
  logs[{EWSudakov_Log_Type::Ls, {}}] = sqr(ls);
  logs[{EWSudakov_Log_Type::lZ, {}}] = ls;
  logs[{EWSudakov_Log_Type::lC, {}}] = ls;
  logs[{EWSudakov_Log_Type::lYuk, {}}] = ls;
  logs[{EWSudakov_Log_Type::lPR, {}}] = ls;
  for (size_t k {0}; k < m_ampls.NumberOfLegs(); ++k)
    for (size_t l {0}; l < k; ++l)
      logs[{EWSudakov_Log_Type::lSSC, {k, l}}] = ls*std::log(std::abs(
          (m_ampls.BaseAmplitude().Leg(k)->Mom()
           + m_ampls.BaseAmplitude().Leg(l)->Mom()).Abs2())/s);

  // calculate
  //   K = (\sum_{i} (1 + 2 Re(\sum_{c} delta_i^c))|M_i|^2)/(\sum_{i} |M_i|^2),
  // where the sum is over the spin configurations; here, we re-write this eq.
  // as
  //   K = 1 + \sum_{c} (2 \sum{i} Re(delta_i^c) |M_i|^2)/(\sum_{i} |M_i|^2)
  //     = 1 + \sum_{c} delta^c
  // which is more convenient since we want to store the c-dependent
  // contributions delta^c with i integrated over (note that c stands for the
  // coeff type).
  auto kfac = 1.0;
  const auto delta_prefactor = m_ewgroupconsts.m_aew/4./M_PI;
  for (const auto& coeffkv : m_coeffs) {
    auto delta_c_num = 0.0;
    for (size_t i {0}; i < m_spinampls[0].size(); ++i) {
      const auto me2 = norm(m_spinampls[0][i]);
      delta_c_num += (coeffkv.second[i] * logs[coeffkv.first]).real() * me2;
    }
    const auto delta_c = 2 * delta_prefactor * delta_c_num / den;
    p_proc->GetMEwgtinfo()->m_ewsudakovkfacdelta[coeffkv.first.first] += delta_c;
    kfac += delta_c;
  }
  EWSudakov_Calculator::m_kfachisto.Insert(kfac);
  return kfac;
}


void EWSudakov_Calculator::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto nspinampls = ampls.size();
  const auto nspins = ampls.GetSpinCombination(0).size();
  assert(nspins == m_ampls.NumberOfLegs());
  m_coeffs.clear();
  for (const auto& key : m_activecoeffs) {
    switch (key) {
      case EWSudakov_Log_Type::Ls:
      case EWSudakov_Log_Type::lZ:
      case EWSudakov_Log_Type::lC:
      case EWSudakov_Log_Type::lYuk:
        m_coeffs[{key, {}}].resize(nspinampls);
        break;
      case EWSudakov_Log_Type::lPR: {
        // NOTE: For the time being we only set this to S. It may however be
        // useful to have a "running" and a "fixed" setting for users.
        m_ewscale2 = m_ampls.MandelstamS();
        m_comixinterface_he.ResetWithEWParameters(
            m_ewgroupconsts.EvolveEWparameters(m_ewscale2));
        m_coeffs[{key, {}}].resize(nspinampls);
        break;
      }
      case EWSudakov_Log_Type::lSSC:
        for (size_t k{ 0 }; k < nspins; ++k)
          for (size_t l{ 0 }; l < k; ++l)
            m_coeffs[{key, {k, l}}].resize(nspinampls);
        break;
    }
  }
  for (size_t i{0}; i < nspinampls; ++i) {
    m_current_me_value = ampls.Get(i);
    if (m_current_me_value == 0.0) {
      continue;
    }
    m_current_spincombination = ampls.GetSpinCombination(i);
    UpdateGolstoneSpincombinationAndMEPrefactor();
    for (const auto& key : m_activecoeffs) {
      switch (key) {
        case EWSudakov_Log_Type::Ls:
          m_coeffs[{key, {}}][i] = LsCoeff();
          break;
        case EWSudakov_Log_Type::lZ:
          m_coeffs[{key, {}}][i] = lsZCoeff();
          break;
        case EWSudakov_Log_Type::lSSC:
          for (size_t k{0}; k < nspins; ++k) {
            for (size_t l{ 0 }; l < k; ++l) {
              // s-channel-related loops will have vanishing log coeffs
              if (k == 1 && l == 0)
                continue;
              if (nspins == 4 && k == 3 && l == 2)
                continue;
              const auto angularkey
                = Coeff_Map_Key{EWSudakov_Log_Type::lSSC, {k, l}};
              m_coeffs[angularkey][i] = lsLogROverSCoeffs({k, l});
            }
          }
          break;
        case EWSudakov_Log_Type::lC:
          m_coeffs[{key, {}}][i] = lsCCoeff();
          break;
        case EWSudakov_Log_Type::lYuk:
          m_coeffs[{key, {}}][i] = lsYukCoeff();
          break;
        case EWSudakov_Log_Type::lPR:
          m_coeffs[{key, {}}][i] = lsPRCoeff();
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
    if (checker.CheckCoeffs(m_coeffs, m_spinampls[0], mandelstam)) {
      THROW(normal_exit, "Finish after checking EW Sudakov coefficients.");
    } else {
      THROW(fatal_error, "EWSudakov coeffs for this process are not equal to"
                         " the results in hep-ph/0010201.");
    }
  }
}

void EWSudakov_Calculator::UpdateGolstoneSpincombinationAndMEPrefactor()
{
  // correct spin index when a longitudinal vector boson is replaced with a
  // scalar using the Goldstone boson equivalence theorem, also calculate
  // (i)^n, the prefactor that accounts for ME^(Z_L^n) -> i^n ME^(chi^n)
  const auto& base_ampl = m_ampls.BaseAmplitude();
  m_current_goldstone_spincombination.clear();
  const auto nspins = m_current_spincombination.size();
  m_current_goldstone_spincombination.reserve(nspins);
  m_current_goldstone_me_prefactor = 1.0;
  for (size_t i{0}; i < nspins; ++i) {
    auto lambda = m_current_spincombination[i];
    if (lambda == 2) {
      lambda = 0;
      if (base_ampl.Leg(i)->Flav().Kfcode() == kf_Z) {
        m_current_goldstone_me_prefactor *= Complex{0.0, 1.0};
      }
    }
    m_current_goldstone_spincombination.push_back(lambda);
  }
}

Coeff_Value EWSudakov_Calculator::LsCoeff()
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  for (size_t i{0}; i < m_current_spincombination.size(); ++i) {
    const Flavour flav{base_ampl.Leg(i)->Flav()};
    const auto diagonal =
        -m_ewgroupconsts.DiagonalCew(flav, m_current_spincombination[i]) / 2.0;
    coeff += diagonal;
    if (flav.Kfcode() == kf_photon || flav.Kfcode() == kf_Z) {
      // special case of neutral transverse gauge bosons, they mix and hence
      // non-diagonal terms appear, cf. e.g. eq. (6.30);
      // assume they are are actually transverse, because we have already
      // replaced longitudinal ones with Goldstone bosons when calling
      // BaseAmplitude() above
      assert(m_current_spincombination[i] != 2);
      const kf_code newkf = (flav.Kfcode() == kf_Z) ? kf_photon : kf_Z;
      const auto prefactor = -m_ewgroupconsts.NondiagonalCew() / 2.0;
      // we can use amplitudes without replacing longitudinal gauge bosons
      // here, because of (i) eq. (3.4) and (ii) Goldstone boson correction
      // terms are diagonal (hence we only need to pass the right flavour above
      // when using DiagonalCew())
      const auto transformed =
          TransformedAmplitudeValue({{i, newkf}}, m_current_spincombination);
      auto amplratio = transformed / m_current_me_value;
      coeff += prefactor * amplratio;
    }
  }
  return coeff;
}

Coeff_Value EWSudakov_Calculator::lsZCoeff()
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  for (size_t i{0}; i < m_current_spincombination.size(); ++i) {
    const Flavour flav{base_ampl.Leg(i)->Flav()};
    const auto couplings =
        m_ewgroupconsts.IZ2(flav, m_current_spincombination[i]);
    for (const auto coupling : couplings) {
      Complex contrib {coupling.second};
      if (coupling.first != flav) {
        // this is a non-diagonal IZ2 term, i.e. we need to take ME ratios into
        // account
        const Leg_Kfcode_Map key{{i, std::abs(coupling.first)}};
        auto transformed =
            TransformedAmplitudeValue(key, m_current_spincombination);
        if (flav.Kfcode() == kf_chi) {
          // we used the Goldstone theorem (3.4) for Z^L -> chi, hence we need
          // to apply a factor of i
          transformed *= Complex{0.0, 1.0};
        }
        const auto amplratio = transformed / m_current_me_value;
        contrib *= amplratio;
      }
      // NOTE: we use 1/m_cw2 = (mZ/mW)^2 for the argument of the logarithm
      coeff += contrib * std::log(1.0 / m_ewgroupconsts.m_cw2);
    }
  }
  return coeff;
}

Coeff_Value
EWSudakov_Calculator::lsLogROverSCoeffs(const Two_Leg_Indizes& indizes)
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  std::vector<Flavour> flavs;
  flavs.reserve(2);
  for (const auto i : indizes) {
    // NOTE: use antiflavours, because the convention is all-incoming in
    // Denner/Pozzorini whereas for Cluster Amplitudes it's all-outgoing
    flavs.push_back(base_ampl.Leg(i)->Flav().Bar());
  }

  // add contribution for each vector boson connecting the leg pairs

  // photon (IA is always diagonal)
  Coeff_Value coeff_A{1.0};
  for (const auto& flav : flavs) {
    coeff_A *= -flav.Charge();
  }
  coeff += 2.0 * coeff_A;

  // Z
  const auto kcouplings =
      m_ewgroupconsts.IZ(flavs[0], m_current_spincombination[indizes[0]]);
  const auto lcouplings =
      m_ewgroupconsts.IZ(flavs[1], m_current_spincombination[indizes[1]]);
  for (const auto kcoupling : kcouplings) {
    for (const auto lcoupling : lcouplings) {
      auto contrib = 2.0 * kcoupling.second * lcoupling.second;
      if (kcoupling.first != flavs[0] || lcoupling.first != flavs[1]) {
        const Leg_Kfcode_Map key{{indizes[0], std::abs(kcoupling.first)},
                                 {indizes[1], std::abs(lcoupling.first)}};
        auto transformed =
            TransformedAmplitudeValue(key, m_current_spincombination);
        for (const auto& flav : flavs) {
          if (flav.Kfcode() == kf_chi) {
            // we used the Goldstone theorem (3.4) for Z^L -> chi, hence we
            // need to apply a factor of i
            transformed *= Complex{0.0, 1.0};
          }
        }
        const auto amplratio = transformed / m_current_me_value;
        contrib *= amplratio;
      }
      coeff += contrib;
    }
  }

  // W
  for (int i{ 0 }; i < 2; ++i) {
    const auto kplus = (i == 0);
    const auto kcouplings = m_ewgroupconsts.Ipm(
        flavs[0], m_current_spincombination[indizes[0]], kplus);
    const auto lcouplings = m_ewgroupconsts.Ipm(
        flavs[1], m_current_spincombination[indizes[1]], !kplus);
    for (const auto kcoupling : kcouplings) {
      for (const auto lcoupling : lcouplings) {
        const Leg_Kfcode_Map key{{indizes[0], std::abs(kcoupling.first)},
                                 {indizes[1], std::abs(lcoupling.first)}};
        auto transformed =
            TransformedAmplitudeValue(key, m_current_spincombination);
        for (const auto& flav : flavs) {
          if (flav.Kfcode() == kf_chi) {
            // we used the Goldstone theorem (3.4) for Z^L -> chi, hence we
            // need to apply a factor of i
            transformed *= Complex{0.0, 1.0};
          }
        }
        const auto amplratio = transformed / m_current_me_value;
        coeff += 2.0*kcoupling.second*lcoupling.second*amplratio;
      }
    }
  }
  return coeff;
}

Coeff_Value EWSudakov_Calculator::lsCCoeff()
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  const auto nspins = m_current_spincombination.size();
  for (size_t i{0}; i < nspins; ++i) {
    const Flavour flav{ base_ampl.Leg(i)->Flav() };
    if (flav.IsFermion()) {
      const auto contrib =
          3.0 / 2.0 *
          m_ewgroupconsts.DiagonalCew(flav, m_current_spincombination[i]);
      coeff += contrib;
    } else if (flav.Kfcode() == kf_Wplus) {
      assert(m_current_spincombination[i] != 2);
      const auto contrib =
          m_ewgroupconsts.DiagonalBew(flav, m_current_spincombination[i]) / 2.0;
      coeff += contrib;
    } else if (flav.Kfcode() == kf_photon || flav.Kfcode() == kf_Z) {
      assert(m_current_spincombination[i] != 2);
      const auto contrib =
          m_ewgroupconsts.DiagonalBew(flav, m_current_spincombination[i]) / 2.0;
      coeff += contrib;
      if (flav.Kfcode() == kf_Z) {
        const auto transformed = TransformedAmplitudeValue(
            {{i, kf_photon}}, m_current_spincombination);
        auto amplratio = transformed / m_current_me_value;
        coeff += m_ewgroupconsts.NondiagonalBew() * amplratio;
      }
    } else if (flav.Kfcode() == kf_chi || flav.Kfcode() == kf_phiplus) {
      const auto contrib = 2.0*m_ewgroupconsts.DiagonalCew(flav, 0);
      coeff += contrib;
    }
  }
  return coeff;
}

Coeff_Value EWSudakov_Calculator::lsYukCoeff()
{
  Coeff_Value coeff{0.0};
  for (size_t i{0}; i < m_current_spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    if (flav.Kfcode() == kf_t || flav.Kfcode() == kf_b) {
      auto contrib = sqr(flav.Mass()/Flavour{kf_Wplus}.Mass());
      if (m_current_spincombination[i] == 0)
        contrib *= 2.0;
      else
        contrib
          += sqr(flav.IsoWeakPartner().Mass()/Flavour{kf_Wplus}.Mass());
      contrib *= -1.0/(8.0*m_ewgroupconsts.m_sw2);
      coeff += contrib;
    } else if (flav.IsVector() && m_current_spincombination[i] == 2) {
      const auto contrib
        = - 3.0/(4.0*m_ewgroupconsts.m_sw2)
        * sqr(Flavour{kf_t}.Mass()/Flavour{kf_Wplus}.Mass());
      coeff += contrib;
    }
  }
  return coeff;
}

Coeff_Value EWSudakov_Calculator::lsPRCoeff()
{
  const auto deno = TransformedAmplitudeValue(
      m_ampls.GoldstoneBosonReplacements(m_current_spincombination),
      m_current_spincombination,
      &m_comixinterface);
  if (deno == 0.0)
    return 0.0;
  const auto he_me = TransformedAmplitudeValue(
      m_ampls.GoldstoneBosonReplacements(m_current_spincombination),
      m_current_spincombination,
      &m_comixinterface_he);
  Coeff_Value coeff = (he_me / deno - 1.0) * 4. * M_PI /
                      log(m_ewscale2 / m_ewgroupconsts.m_mw2) /
                      m_ewgroupconsts.m_aew;
  return coeff;
}

Complex EWSudakov_Calculator::TransformedAmplitudeValue(
    const Leg_Kfcode_Map& legs,
    const std::vector<int>& spincombination,
    const Comix_Interface* interface)
{
  auto& transformedampl = m_ampls.SU2TransformedAmplitude(legs);
  /// TODO: Make the following a bit prettier. At the moment
  /// this is simply a flag to catch that the momentum
  /// stretcher has failed. May want to have a enum
  if (transformedampl.Flag() & (1 << 4)) {
    ++m_numonshellwarning;
    return 0.0;
  }

  // correct for Goldstone boson replacements
  std::vector<int> guildedspincombination;
  size_t i {0};
  for (int i {0}; i < spincombination.size(); ++i) {
    auto pol = spincombination[i];
    if (pol == 2) {
      auto it = legs.find(i);
      if (it != legs.end()) {
        if (it->second == kf_chi || it->second == kf_phiplus ||
            it->second == kf_h0) {
          pol = 0;
        }
      }
    }
    guildedspincombination.push_back(pol);
  }

  return (interface ? *interface : m_comixinterface)
      .GetSpinAmplitude(transformedampl, guildedspincombination);
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
