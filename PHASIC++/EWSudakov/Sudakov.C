#include "PHASIC++/EWSudakov/Sudakov.H"

#include "PHASIC++/EWSudakov/Comix_Interface.H"
#include "PHASIC++/EWSudakov/Coefficient_Checker.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/KF_Table.H"

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
    EWSudakov_Log_Type::lYuk,
    EWSudakov_Log_Type::lPR
  },
  m_ampls{ p_proc, m_activecoeffs },
  m_comixinterface{ p_proc, m_ampls },
  m_mw2{ sqr(s_kftable[kf_Wplus]->m_mass) },
  m_runaqed{ 1./137.03599976 }
{
  auto& s = Settings::GetMainSettings();
  m_check = s["CHECK_EWSUDAKOV"].SetDefault(false).Get<bool>();
  m_threshold = s["EWSUDAKOV_THRESHOLD"].SetDefault(5.0).Get<double>();
}

double Sudakov::KFactor(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC("");
  m_ampls.UpdateMomenta(mom);
  if (!IsInHighEnergyLimit())
    return 1.0;
  ClearSpinAmplitudes();
  FillBaseSpinAmplitudes();
  CalculateSpinAmplitudeCoeffs();
  return KFactor();
}

bool Sudakov::IsInHighEnergyLimit()
{
  DEBUG_FUNC("");
  static const auto threshold = sqr(m_threshold) * m_mw2;

  const auto& base_ampl = m_ampls.BaseAmplitude();
  for (size_t i {0}; i < base_ampl.Legs().size(); ++i) {
    for (size_t j {i + 1}; j <  base_ampl.Legs().size(); ++j) {
      const auto sij
        = std::abs((base_ampl.Leg(i)->Mom() + base_ampl.Leg(j)->Mom()).Abs2());
      if(sij < threshold)
        return false;
    }
  }
  return true;
}

void Sudakov::ClearSpinAmplitudes()
{
  m_spinampls.clear();
  m_transformedspinampls.clear();
}

void Sudakov::FillBaseSpinAmplitudes()
{
  m_comixinterface.FillSpinAmplitudes(m_spinampls, m_ampls.BaseAmplitude());
}

double Sudakov::KFactor()
{
  auto den = m_spinampls[0].SumSquare();
  if (den == 0.0)
    return 1.0;
  const auto s = std::abs(m_ampls.MandelstamS());
  const auto ls = std::log(s/m_mw2);

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

  // calculate K = (\sum_{i} (1 + 2 Re(delta))|M_i|^2) / (\sum_{i} |M_i|^2),
  // where the sum is over the spin configurations
  auto num = 0.0;
  for (size_t i {0}; i < m_spinampls[0].size(); ++i) {
    static const auto delta_prefactor = m_runaqed.AqedThomson()/4./M_PI;
    auto delta = 0.0;
    for (const auto& coeffkv : m_coeffs)
      delta += (coeffkv.second[i].first * logs[coeffkv.first]).real();
    num += (1.0 + 2.0*delta_prefactor*delta) * norm(m_spinampls[0][i]);
  }

  return num/den;
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
      case EWSudakov_Log_Type::lPR:
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
        case EWSudakov_Log_Type::lPR:
          m_coeffs[{key, {}}][i] = lsPRCoeff(value, spincombination);
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

Coeff_Value Sudakov::LsCoeff(Complex amplvalue,
                             std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i{ 0 }; i < spincombination.size(); ++i) {
    const Flavour flav{m_ampls.BaseAmplitude(spincombination).Leg(i)->Flav()};
    const auto diagonal
      = -m_ewgroupconsts.DiagonalCew(flav, spincombination[i]) / 2.0;
    coeff.first += diagonal;
    coeff.second += diagonal;
    if (flav.IsVector() && flav.Charge() == 0 && spincombination[i] != 2) {
      const kf_code newkf = (flav.Kfcode() == kf_Z) ? kf_photon : kf_Z;
      // special case of neutral gauge bosons, they mix and hence non-diagonal
      // terms appear, cf. e.g. eq. (6.30)
      const auto prefactor = -m_ewgroupconsts.NondiagonalCew() / 2.0;
      const auto transformed
        = TransformedAmplitudeValue({{i, newkf}}, spincombination);
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
  // NOTE: use antiflavours, because the convention is all-incoming in
  // Denner/Pozzorini whereas for Cluster Amplitudes it's all-outgoing
  auto kflav = m_ampls.BaseAmplitude().Leg(k)->Flav().Bar();
  auto lflav = m_ampls.BaseAmplitude().Leg(l)->Flav().Bar();

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
        const Leg_Kfcode_Map key{{k, kcoupling.first}, {l, lcoupling.first}};
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
        const auto transformed
          = TransformedAmplitudeValue(key, goldstonespincombination);
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
    if (flav.IsFermion()) {
      const auto contrib
        = 3.0/2.0 * m_ewgroupconsts.DiagonalCew(flav, spincombination[i]);
      coeff.first += contrib;
      coeff.second += contrib;
    } else if (flav.Kfcode() == kf_Wplus && spincombination[i] != 2) {
      const auto contrib
        = m_ewgroupconsts.DiagonalBew(flav, spincombination[i]) / 2.0;
      coeff.first += contrib;
      coeff.second += contrib;
    } else if (flav.IsVector()
               && flav.Charge() == 0
               && spincombination[i] != 2) {
      const auto contrib
        = m_ewgroupconsts.DiagonalBew(flav, spincombination[i]) / 2.0;
      coeff.first += contrib;
      coeff.second += contrib;
      if (flav.Kfcode() == kf_Z) {
        const auto transformed
          = TransformedAmplitudeValue({{i, kf_photon}}, spincombination);
        const auto base = amplvalue;
        assert(base != 0.0);  // guaranteed by CalculateSpinAmplitudeCoeffs
        auto amplratio = transformed/base;
        coeff.first += m_ewgroupconsts.NondiagonalBew() * amplratio;
        coeff.second -= m_ewgroupconsts.NondiagonalBew() * amplratio;
      }
    } else if (flav.IsVector() && spincombination[i] == 2) {
      const auto contrib = 2.0*m_ewgroupconsts.DiagonalCew(flav, 2);
      coeff.first += contrib;
      coeff.second += contrib;
    }
  }
  return coeff;
}

Coeff_Value Sudakov::lsYukCoeff(Complex amplvalue,
                                std::vector<int> spincombination)
{
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  for (size_t i {0}; i < spincombination.size(); ++i) {
    const Flavour flav{ m_ampls.BaseAmplitude().Leg(i)->Flav() };
    if (flav.Kfcode() == kf_t || flav.Kfcode() == kf_b) {
      auto contrib = sqr(flav.Mass()/Flavour{kf_Wplus}.Mass());
      if (spincombination[i] == 0)
        contrib *= 2.0;
      else
        contrib
          += sqr(flav.IsoWeakPartner().Mass()/Flavour{kf_Wplus}.Mass());
      contrib *= -1.0/(8.0*m_ewgroupconsts.m_sw2);
      coeff.first += contrib;
      coeff.second += contrib;
    } else if (flav.IsVector() && spincombination[i] == 2) {
      const auto contrib
        = - 3.0/(4.0*m_ewgroupconsts.m_sw2)
        * sqr(Flavour{kf_t}.Mass()/Flavour{kf_Wplus}.Mass());
      coeff.first += contrib;
      coeff.second += contrib;
    }
  }
  return coeff;
}

Coeff_Value Sudakov::lsPRCoeff(Complex amplvalue,
                               std::vector<int> spincombination)
{
  // TODO: implement dynamic calculation, this placeholder just returns the
  // coefficients given in the Denner/Pozzorini reference for ee->mumu
  auto coeff = std::make_pair(Complex{ 0.0 }, Complex{ 0.0 });
  if (spincombination[0] == 0
      && spincombination[1] == 0
      && spincombination[2] == 0
      && spincombination[3] == 0)
    coeff.first = coeff.second = 8.80;
  else if (spincombination[0] == 0
      && spincombination[1] == 0
      && spincombination[2] == 1
      && spincombination[3] == 1)
    coeff.first = coeff.second = 8.80;
  else if (spincombination[0] == 1
      && spincombination[1] == 1
      && spincombination[2] == 0
      && spincombination[3] == 0)
    coeff.first = coeff.second = 8.80;
  else if (spincombination[0] == 1
      && spincombination[1] == 1
      && spincombination[2] == 1
      && spincombination[3] == 1)
    coeff.first = coeff.second = -9.03;
  return coeff;
}

Complex
Sudakov::TransformedAmplitudeValue(const Leg_Kfcode_Map& legs,
                                   const std::vector<int>& spincombination)
{
  auto amplit = m_transformedspinampls.find(legs);
  if (amplit == m_transformedspinampls.end()) {
    auto& transformedampl = m_ampls.SU2TransformedAmplitude(legs);
    m_comixinterface.FillSpinAmplitudes(m_transformedspinampls[legs],
                                        transformedampl);
    amplit = m_transformedspinampls.find(legs);
  }
  auto& legpermutation = m_ampls.LegPermutation(legs);
  std::vector<int> transformedspincombination;
  for (const auto& idx : legpermutation)
    transformedspincombination.push_back(spincombination[idx]);
  return amplit->second[0].Get(transformedspincombination);
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
