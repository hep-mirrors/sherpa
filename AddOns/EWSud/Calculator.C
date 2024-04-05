#include "AddOns/EWSud/Calculator.H"

#include "Coefficient_Checker.H"
#include "KFactor_Checker.H"
#include "PHASIC++/Main/Process_Integrator.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <cassert>

using namespace PHASIC;
using namespace ATOOLS;
using namespace EWSud;

Histogram Calculator::m_kfachisto(0, -5.0, 5.0, 50);
size_t Calculator::m_numonshellwarning {0};

Calculator::Calculator(Process_Base* proc):
  p_proc{ proc },
  m_activelogtypes{
    EWSudakov_Log_Type::Ls,
    EWSudakov_Log_Type::lZ,
    EWSudakov_Log_Type::lSSC,
    EWSudakov_Log_Type::lC,
    EWSudakov_Log_Type::lYuk,
    EWSudakov_Log_Type::lPR,
    EWSudakov_Log_Type::lI
  },
  m_ampls{ p_proc, m_activelogtypes },
  m_comixinterface{ p_proc, m_ampls },
  m_comixinterface_he{ p_proc, m_ampls }
{

  Scoped_Settings s = Settings::GetMainSettings()["EWSUD"];
  m_checkcoeff = s["CHECK"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_CHECK"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: CHECK");
  }

  m_checkkfac = s["CHECK_KFACTOR"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_CHECK_KFACTOR"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: CHECK_KFACTOR");
  }

  m_checklogfile =
      s["CHECK_LOG_FILE"].SetDefault("").Get<std::string>();
  if(Settings::GetMainSettings()["CHECK_EWSUDAKOV_LOG_FILE"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: CHECK_EWSUDAKOV_LOG_FILE");
  }

  m_threshold = s["THRESHOLD"].SetDefault(5.0).Get<double>();
  if(Settings::GetMainSettings()["EWSUDAKOV_THRESHOLD"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD:EWSUDAKOV_THRESHOLD");
  }

  m_checkinvariantratios = s["CHECKINVARIANTRATIOS"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_CHECKINVARIANTRATIOS"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: EWSUDAKOV_CHECKINVARIANTRATIOS");
  }

  s.DeclareVectorSettingsWithEmptyDefault({"COEFF_REMOVED_LIST"});
  if(Settings::GetMainSettings()["EWSUDAKOV_COEFF_REMOVED_LIST"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: EWSUDAKOV_COEFF_REMOVED_LIST");
  }

  const auto disabled_log_list =
      s["COEFF_REMOVED_LIST"].GetVector<std::string>();
  if(Settings::GetMainSettings()["EWSUDAKOV_COEFF_REMOVED_LIST"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD:EWSUDAKOV_COEFF_REMOVED_LIST");
  }

  for (const auto& l : disabled_log_list) {
    m_activelogtypes.erase(EWSudakovLogTypeFromString(l));
  }
  m_c_coeff_ignores_vector_bosons =
      s["C_COEFF_IGNORES_VECTOR_BOSONS"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_C_COEFF_IGNORES_VECTOR_BOSONS"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: C_COEFF_IGNORES_VECTOR_BOSONS");
  }

  SetHighEnergyScheme(s["HIGH_ENERGY_SCHEME"].SetDefault("Default").Get<std::string>());
  if(Settings::GetMainSettings()["EWSUDAKOV_HIGH_ENERGY_SCHEME"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: HIGH_ENERGY_SCHEME");
  }

  m_includesubleading = s["INCLUDE_SUBLEADING"].SetDefault(false).Get<bool>();
  if(Settings::GetMainSettings()["EWSUDAKOV_INCLUDE_SUBLEADING"].IsSetExplicitly()){
    THROW(fatal_error, "Avoid Using old syntax, prefer the new EWSUD: INCLUDE_SUBLEADING");
  }

}

Calculator::~Calculator()
{
  static bool did_output{false};
  if (!did_output) {
    Calculator::m_kfachisto.MPISync();
    Calculator::m_kfachisto.Finalize();
    MyStrStream s;
    s << "kfacs_" << m_threshold;
    Calculator::m_kfachisto.Output(s.str());
    msg_Error() << "Set " << m_numonshellwarning
                << " amplitudes to 0.0, because there was not enough energy to "
                   "fulfil on-shell conditions\n";
    did_output = true;
  }
}

void Calculator::SetHighEnergyScheme(const std::string& hescheme)
{
  if(hescheme == "Default"){
    m_helimitscheme = HighEnergySchemes::_default;
  } else if (hescheme == "Tolerant"){
    m_helimitscheme = HighEnergySchemes::_tolerant;
  } else if (hescheme == "Cluster_Dumb"){
    m_helimitscheme = HighEnergySchemes::_cluster_dumb;
  } else if (hescheme == "Cluster_L1"){
    m_helimitscheme = HighEnergySchemes::_cluster_l1;
  } else {
    THROW(not_implemented, ("Option " + hescheme
      + " is not implemented, valid options are: Default, Tolerant, Cluster_Dumb and Cluster_L1"));
  }
};

EWSudakov_Log_Corrections_Map
Calculator::CorrectionsMap(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);

#if false
  if (subcalculator) {
    Vec4D_Vector cmom = mom;
    cmom[2] = cmom[2] + cmom[3];
    cmom.erase(cmom.begin() + 3);
    return subcalculator->CorrectionsMap(cmom);
  } else {
#endif
    if (!IsInHighEnergyLimit())
      return {};
#if false
  }
#endif

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
  return CorrectionsMap();
}

bool Calculator::IsInHighEnergyLimit()
{
  DEBUG_FUNC("");
  /// In all schemes but the default one, do the checks when actually
  /// calculating the logarithms
  if(m_helimitscheme != HighEnergySchemes::_default) return true;
  static const auto threshold = sqr(m_threshold) * m_ewgroupconsts.m_mw2;

  const auto s = std::abs(m_ampls.MandelstamS());

  const auto& base_ampl = m_ampls.BaseAmplitude();
  for (size_t i {0}; i < base_ampl.Legs().size(); ++i) {
    for (size_t j {i + 1}; j <  base_ampl.Legs().size(); ++j) {
      const auto sij
        = std::abs((base_ampl.Mom(i) + base_ampl.Mom(j)).Abs2());
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

void Calculator::ClearSpinAmplitudes()
{
  m_spinampls.clear();
  m_comixinterface.ResetSpinAmplitudeCache();
  m_comixinterface_he.ResetSpinAmplitudeCache();
}

void Calculator::FillBaseSpinAmplitudes()
{
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.BaseAmplitude());
}

EWSudakov_Log_Corrections_Map Calculator::CorrectionsMap()
{
  auto den = m_spinampls[0].SumSquare();
  if (den == 0.0) {
    return {};
  }
  const auto s = std::abs(m_ampls.MandelstamS());
  const auto ls = std::log(s/m_ewgroupconsts.m_mw2);

  // pre-calculate the logarithms we need below
  std::map<Coeff_Map_Key, double> logs;
  logs[{EWSudakov_Log_Type::Ls, {}}] = sqr(ls);
  logs[{EWSudakov_Log_Type::lZ, {}}] = ls;
  logs[{EWSudakov_Log_Type::lC, {}}] = ls;
  logs[{EWSudakov_Log_Type::lYuk, {}}] = ls;
  logs[{EWSudakov_Log_Type::lPR, {}}] = ls;
  const auto& base_ampl = m_ampls.BaseAmplitude();
  for (size_t k{0}; k < m_ampls.NumberOfLegs(); ++k) {
    for (size_t l{0}; l < k; ++l) {
      /// If the scheme is not the default, I still need to perform checks
      /// on the invariants
      if (m_helimitscheme == HighEnergySchemes::_tolerant ) {
        const double rkl       = std::abs((base_ampl.Mom(k) + base_ampl.Mom(l)).Abs2());
        const double threshold = sqr(m_threshold) * m_ewgroupconsts.m_mw2;
        logs[{EWSudakov_Log_Type::lSSC, {k, l}}] =
          (rkl < threshold) ? 0.0 : ls * std::log(rkl / s);
      } else { /// for the moment I don't have the other schemes
        const double rkl =
            std::abs((base_ampl.Mom(k) + base_ampl.Mom(l)).Abs2());
        double log = ls * std::log(rkl / s);
        if (m_includesubleading)
          log += sqr(std::log(rkl / s)) / 2.0; // last term in eq. (3.3)
        logs[{EWSudakov_Log_Type::lSSC, {k, l}}] = log;
      }
    }
  }

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
  EWSudakov_Log_Corrections_Map kfacs;
  const auto delta_prefactor = m_ewgroupconsts.delta_prefactor;
  for (const auto& coeffkv : m_coeffs) {
    auto delta_c_num = 0.0;
    for (size_t i {0}; i < m_spinampls[0].size(); ++i) {
      const auto me2 = norm(m_spinampls[0][i]);
      delta_c_num += (coeffkv.second[i] * logs[coeffkv.first]).real() * me2;
    }
    const auto delta_c = 2 * delta_prefactor * delta_c_num / den;
    kfacs[coeffkv.first.first] += delta_c;
    kfac += delta_c;
  }
  Calculator::m_kfachisto.Insert(kfac);
  if (m_checkkfac) {
    KFactor_Checker checker(p_proc->Name());
    checker.SetLogFileName(m_checklogfile);
    Mandelstam_Variables mandelstam {
      m_ampls.MandelstamS(),
      m_ampls.MandelstamT(),
      m_ampls.MandelstamU() };
    if (checker.CheckKFactor(
            kfacs.KFactor(), mandelstam, m_ewgroupconsts)) {
      THROW(normal_exit, "Finish after checking EWSud K factor.");
    } else {
      THROW(fatal_error, "EWSud K factor for this process is not equal to"
                         " the results in hep-ph/0408308.");
    }
  }
  return kfacs;
}


void Calculator::CalculateSpinAmplitudeCoeffs()
{
  const auto& ampls = m_spinampls[0];
  const auto nspinampls = ampls.size();
  const auto nspins = ampls.GetSpinCombination(0).size();
  assert(nspins == m_ampls.NumberOfLegs());
  m_coeffs.clear();
  for (const auto& key : m_activelogtypes) {
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
      case EWSudakov_Log_Type::lI:
        // This is calculated *before* potential clusterings, i.e. in the
        // Clustered_Calculator
        break;
    }
  }
  for (size_t i{0}; i < nspinampls; ++i) {
    m_current_me_value = ampls.Get(i);
    if (m_current_me_value == 0.0) {
      continue;
    }
    m_current_spincombination = ampls.GetSpinCombination(i);
    for (const auto& key : m_activelogtypes) {
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
        case EWSudakov_Log_Type::lI:
          // This is calculated *before* potential clusterings, i.e. in the
          // Clustered_Calculator
          break;
      }
    }
  }
  if (m_checkcoeff) {
    Coefficient_Checker checker(p_proc->Name(), m_activelogtypes);
    checker.SetLogFileName(m_checklogfile);
    Mandelstam_Variables mandelstam {
      m_ampls.MandelstamS(),
      m_ampls.MandelstamT(),
      m_ampls.MandelstamU() };
    if (checker.CheckCoeffs(
            m_coeffs, m_spinampls[0], mandelstam, m_ewgroupconsts)) {
      THROW(normal_exit, "Finish after checking EWSud coefficients.");
    } else {
      THROW(fatal_error, "EWSud coeffs for this process are not equal to"
                         " the results in hep-ph/0010201.");
    }
  }
}

Complex Calculator::GBETConversionFactor(Leg_Kfcode_Map legs)
{
  // calculate the Goldstone boson equivalence theorem factor for going from
  // the physical phase to the Goldstone phase, for the current amplitude after
  // applying any leg changes encoded in `legs`

  // the factors are taken from arXiv:hep-ph/0201077 (Pozzorini's thesis),
  // I (EB) was not able to use it to calculate ME ratios in the physical phase
  // for the longitudinal Goldstone bosons in the various calculator functions
  // as e.g. lsLogROverSCoeffs. Even when combined with the conversion factors
  // given below (and various variants of that), the ME ratios times the
  // conversion factors were not equal to the ME ratios in the Goldstone phase,
  // which seem to give the correct result in all test cases; I leave the
  // implementation for future reference, but throw to prevent usage
  THROW(fatal_error, "Called GBETConversionFactor(), which is likely wrong.");

  // add longitudinal boson -> Goldstone boson replacements and get the ampl
  auto tmp = m_ampls.GoldstoneBosonReplacements(m_current_spincombination);
  legs.insert(std::begin(tmp), std::end(tmp));
  if (legs.empty())
    return 1.0;
  auto& ampl = m_ampls.SU2TransformedAmplitude(legs);

  const auto nspins = m_current_spincombination.size();
  Complex factor = 1.0;
  for (size_t i{0}; i < nspins; ++i) {
    if (ampl.Leg(i)->Flav().Kfcode() == kf_chi) {
      factor *= Complex{0.0, 1.0};
    } else if (ampl.Leg(i)->Flav().Kfcode() == kf_phiplus && ampl.Leg(i)->Flav().Bar().IsAnti()) {
      factor *= Complex {-1.0, 0.0};
    }
  }
  return factor;
}

Coeff_Value Calculator::LsCoeff()
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
      const auto transformed =
          TransformedAmplitudeValue({{i, newkf}}, m_current_spincombination);
      const auto amplratio = transformed / m_current_me_value;
      coeff += prefactor * amplratio;
    }
  }
  return coeff;
}

Coeff_Value Calculator::lsZCoeff()
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  for (size_t i{0}; i < m_current_spincombination.size(); ++i) {
    const Flavour flav{base_ampl.Leg(i)->Flav().Bar()};
    const auto IZ2 = m_ewgroupconsts.IZ2(flav, m_current_spincombination[i]);
    // NOTE: we use 1/m_cw2 = (mZ/mW)^2 for the argument of the logarithm
    coeff += IZ2 * std::log(1.0 / m_ewgroupconsts.m_cw2);
  }
  return coeff;
}

Coeff_Value
Calculator::lsLogROverSCoeffs(const Two_Leg_Indizes& indizes)
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
    coeff_A *= flav.Charge();
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
      if (Flavour(kcoupling.first) != flavs[0] || Flavour(lcoupling.first) != flavs[1]) {
        Leg_Kfcode_Map key {{indizes[0], std::abs(kcoupling.first)},
                            {indizes[1], std::abs(lcoupling.first)}};
        const auto goldstone_legs = m_ampls.GoldstoneBosonReplacements(m_current_spincombination);
        key.insert(goldstone_legs.begin(), goldstone_legs.end());
        const auto transformed =
            TransformedAmplitudeValue(key, m_current_spincombination);
        const auto deno = TransformedAmplitudeValue(
            m_ampls.GoldstoneBosonReplacements(m_current_spincombination),
            m_current_spincombination,
            &m_comixinterface);
        // NOTE: deno can be zero due to a process not being found, e.g.
        // g dbar -> Z dbar ~> g dbar -> \chi dbar, but YUKAWA set to zero
        const auto amplratio = (deno == 0.0) ? 0.0 : transformed / deno;
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
        Leg_Kfcode_Map key {{indizes[0], std::abs(kcoupling.first)},
                            {indizes[1], std::abs(lcoupling.first)}};
        const auto goldstone_legs = m_ampls.GoldstoneBosonReplacements(m_current_spincombination);
        key.insert(goldstone_legs.begin(), goldstone_legs.end());
        const auto transformed =
            TransformedAmplitudeValue(key, m_current_spincombination);
        const auto deno = TransformedAmplitudeValue(
            m_ampls.GoldstoneBosonReplacements(m_current_spincombination),
            m_current_spincombination,
            &m_comixinterface);
        // NOTE: deno can be zero, cf. comment above for Z loop terms
        const auto amplratio = (deno == 0.0) ? 0.0 : transformed / deno;
        coeff += 2.0*kcoupling.second*lcoupling.second*amplratio;
      }
    }
  }
  return coeff;
}

Coeff_Value Calculator::lsCCoeff()
{
  Coeff_Value coeff{0.0};
  const auto& base_ampl = m_ampls.BaseAmplitude(m_current_spincombination);
  const auto nspins = m_current_spincombination.size();
  for (size_t i{0}; i < nspins; ++i) {
    const Flavour flav{ base_ampl.Leg(i)->Flav().Bar() };
    if (flav.IsFermion()) {
      const auto contrib =
          3.0 / 2.0 *
          m_ewgroupconsts.DiagonalCew(flav, m_current_spincombination[i]);
      coeff += contrib;
    } else if (flav.Kfcode() == kf_Wplus) {
      assert(m_current_spincombination[i] != 2);
      if (m_c_coeff_ignores_vector_bosons) continue;
      const auto contrib =
          m_ewgroupconsts.DiagonalBew(flav, m_current_spincombination[i]) / 2.0;
      coeff += contrib;
    } else if (flav.Kfcode() == kf_photon || flav.Kfcode() == kf_Z) {
      assert(m_current_spincombination[i] != 2);
      if (m_c_coeff_ignores_vector_bosons) continue;
      const auto contrib =
          m_ewgroupconsts.DiagonalBew(flav, m_current_spincombination[i]) / 2.0;
      coeff += contrib;
      if (flav.Kfcode() == kf_Z) {
        const auto transformed =
            TransformedAmplitudeValue({{i, kf_photon}},
                                      m_current_spincombination);
        const auto amplratio = transformed / m_current_me_value;
        coeff += m_ewgroupconsts.NondiagonalBew() * amplratio;
      }
    } else if (flav.Kfcode() == kf_chi || flav.Kfcode() == kf_phiplus) {
      const auto contrib = 2.0*m_ewgroupconsts.DiagonalCew(flav, 0);
      coeff += contrib;
    }
  }
  return coeff;
}

Coeff_Value Calculator::lsYukCoeff()
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

Coeff_Value Calculator::lsPRCoeff()
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

Complex Calculator::TransformedAmplitudeValue(
    const Leg_Kfcode_Map& legs,
    const std::vector<int>& spincombination,
    const Comix_Interface* interface)
{
  auto& transformedampl = m_ampls.SU2TransformedAmplitude(legs);
  if (transformedampl.Flag() == m_ampls.StretcherFailFlag) {
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


namespace EWSud {

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
