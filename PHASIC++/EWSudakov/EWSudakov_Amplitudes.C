#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/EWSudakov/EWGroupConstants.H"
#include "PHASIC++/Process/Process_Base.H"

#include <numeric>
#include <cassert>

using namespace PHASIC;
using namespace ATOOLS;

const EWSudakov_Amplitudes::Cluster_Ampl_Key EWSudakov_Amplitudes::s_baseamplkey
= Leg_Kfcode_Map{};

EWSudakov_Amplitudes::EWSudakov_Amplitudes(
    Process_Base* proc, const std::set<EWSudakov_Log_Type>& activecoeffs):
  ampls{ CreateAmplitudes(proc, activecoeffs) },
  permutations{ CreatePermutations(ampls) }
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::BaseAmplitude() noexcept
{
  return SU2TransformedAmplitude(s_baseamplkey);
}

Cluster_Amplitude&
EWSudakov_Amplitudes::BaseAmplitude(std::vector<int> spincombination)
{
  Leg_Kfcode_Map leg_set;
  for (int i{ 0 }; i < BaseAmplitude().Legs().size(); ++i) {
    if (spincombination[i] == 2) {
      const auto flav = BaseAmplitude().Leg(i)->Flav();
      if (flav.Kfcode() == kf_Z || flav.Kfcode() == kf_Wplus) {
        leg_set.emplace(i, flav.GoldstoneBosonPartner().Kfcode());
      }
    }
  }
  return SU2TransformedAmplitude(leg_set);
}

Cluster_Amplitude& EWSudakov_Amplitudes::SU2TransformedAmplitude(
  const Leg_Kfcode_Map& legs)
{
  const auto it = ampls.find(legs);
  if (it == ampls.end()) {
    MyStrStream s;
    s << "SU(2)-transformed amplitude not found:\n" << legs;
    THROW(fatal_error, s.str());
  }
  return *(it->second);
}

std::vector<size_t>& EWSudakov_Amplitudes::LegPermutation(const Leg_Kfcode_Map& legs)
{
  const auto it = permutations.find(legs);
  if (it == permutations.end())
    THROW(fatal_error, "Permutation not found");
  return it->second;
}

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC(mom);
  for (int i{ 0 }; i < BaseAmplitude().Legs().size(); ++i) {
    BaseAmplitude().Leg(i)->SetMom(i < BaseAmplitude().NIn() ? -mom[i] : mom[i]);
  }
  // TODO: generalise the momentum rescaling below to handle 2->n
  assert(BaseAmplitude().NIn() == 2 && BaseAmplitude().Legs().size() == 4);
  // TODO: also, reuse existing implementation elsewhere if existing
  // TODO: drop the implicit assumption that all initial-state particles are
  // massless (?), e.g. there might be top-quarks in the IS after SU(2)
  // transforming the amplitude
  const auto Etot2 = MandelstamS();

  // boost to the rest system if necessary
  bool should_boost {false};
  Poincare rest;
  Vec4D cms = Vec4D {0.0, 0.0, 0.0, 0.0};
  for (int i {0}; i < 2; i++) {
    cms += BaseAmplitude().Leg(i)->Mom();
  }
  if (Vec3D(cms).Abs() > 1.e-6) {
    should_boost = true;
    rest = Poincare(cms);
  }

  for (auto& ampl : ampls) {
    if (ampl.first == s_baseamplkey)
      continue;
    const auto& permutation = LegPermutation(ampl.first);

    // first modify incoming momenta, then outgoing momenta
    for (int i {0}; i < 4; i += 2) {
      // calc momenta such that all particles are on their mass shell;
      // to guarantee this we modify the energy and the absolute value of the
      // momenta, but not the directions of the momenta
      Vec4D mom0 {BaseAmplitude().Leg(permutation[i])->Mom()};
      if (should_boost) {
        rest.Boost(mom0);
      }
      const Vec3D normed_mom0 {Vec3D {mom0} / Vec3D {mom0}.Abs()};
      const auto m02 =
          sqr(ATOOLS::Abs<double>(ampl.second->Leg(i)->Flav().Mass()));
      const auto m12 =
          sqr(ATOOLS::Abs<double>(ampl.second->Leg(i + 1)->Flav().Mass()));
      // calc the fraction of energy that will be assigned to the first leg
      const auto x = (Etot2 + m02 - m12) / (2 * Etot2);
      // calc the absolute momentum assigned to the first leg
      const auto p =
          sqrt((Etot2 * Etot2 - 2 * Etot2 * (m02 + m12) + sqr(m02 - m12)) /
               (4 * Etot2));

      mom0 = Vec4D {x * sqrt(Etot2), p * normed_mom0};
      Vec4D mom1 {(1 - x) * sqrt(Etot2), -p * normed_mom0};

      if (should_boost) {
        rest.BoostBack(mom0);
        rest.BoostBack(mom1);
      }

      if (i == 0) {
        mom0 = -mom0;
        mom1 = -mom1;
      }

      ampl.second->Leg(i)->SetMom(mom0);
      ampl.second->Leg(i + 1)->SetMom(mom1);
    }
  }
}

double EWSudakov_Amplitudes::MandelstamS()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Leg(0)->Mom() + ampl.Leg(1)->Mom()).Abs2();
}

double EWSudakov_Amplitudes::MandelstamT()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Leg(0)->Mom() - ampl.Leg(2)->Mom()).Abs2();
}

double EWSudakov_Amplitudes::MandelstamU()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Leg(0)->Mom() - ampl.Leg(3)->Mom()).Abs2();
}

EWSudakov_Amplitudes::Cluster_Amplitude_UPM
EWSudakov_Amplitudes::CreateAmplitudes(
    Process_Base* proc,
    const std::set<EWSudakov_Log_Type>& activecoeffs) const
{
  Cluster_Amplitude_UPM ampls;

  const EWGroupConstants ewgroupconsts;

  // create unmodified amplitude
  const auto& baseampl =
      ampls.insert(std::make_pair(s_baseamplkey, CreateAmplitude(proc)))
          .first->second;
  const auto nlegs = baseampl->Legs().size();

  // iterate over permutations of Z -> \chi and W -> phi Goldstone boson
  // replacements

  // store W and Z indizes
  std::vector<size_t> bosonindexes;
  bosonindexes.reserve(nlegs);
  for (size_t k{0}; k < nlegs; ++k) {
    const auto flav = baseampl->Leg(k)->Flav();
    if (flav.Kfcode() == kf_Z || flav.Kfcode() == kf_Wplus) {
      bosonindexes.push_back(k);
    }
  }

  // permute over replacing / not replacing each boson index
  const size_t first_invalid_permutation{
      static_cast<size_t>(1 << bosonindexes.size())};
  for (size_t permutation{0}; permutation != first_invalid_permutation;
       ++permutation) {
    auto* current_ampl{&baseampl};
    if (permutation != 0) {
      Leg_Kfcode_Map leg_set;
      Leg_Kfcode_Map_Signed leg_set_signed;
      for (size_t k{0}; k < bosonindexes.size(); ++k) {
        if (permutation & (1 << k)) {
          const auto bosonindex = bosonindexes[k];
          const auto& flav = (*current_ampl)->Leg(bosonindex)->Flav();
          const auto flavcode = (long int)flav.GoldstoneBosonPartner();
          leg_set.emplace(bosonindex, std::abs(flavcode));
          leg_set_signed.emplace(bosonindex, flavcode);
        }
      }
      auto it = ampls.find(leg_set);
      if (it == ampls.end()) {
        auto ampl = std::make_pair(
            leg_set, CreateSU2TransformedAmplitude((*current_ampl), leg_set_signed));
        it = ampls.insert(std::move(ampl)).first;
      }
      current_ampl = &it->second;
    }

    // create ampls needed for Z/photon mixing in Ls coefficients (induced by
    // non-diagonal elements of C^ew)
    if (activecoeffs.find(EWSudakov_Log_Type::Ls) != activecoeffs.end()) {
      for (size_t i{0}; i < nlegs; ++i) {
        const auto flav = (*current_ampl)->Leg(i)->Flav();
        int newkf{kf_none};
        if (flav.IsPhoton())
          newkf = kf_Z;
        else if (flav.Kfcode() == kf_Z)
          newkf = kf_photon;
        if (newkf != kf_none) {
          auto ampl = std::make_pair(
              Cluster_Ampl_Key{{i, newkf}},
              CreateSU2TransformedAmplitude(
                  (*current_ampl), {{i, static_cast<long int>(newkf)}}));
          ampls.insert(std::move(ampl));
        }
      }
    }

    // create ampls needed for terms with a single I^Z (either squared or
    // non-squared, i.e. it does not matter if we use IZ or IZ2 below because
    // we only need the corresponding flavour replacement)
    if (activecoeffs.find(EWSudakov_Log_Type::lZ) != activecoeffs.end()) {
      for (size_t i{0}; i < nlegs; ++i) {
        const auto flav = (*current_ampl)->Leg(i)->Flav();
        const auto couplings = ewgroupconsts.IZ2(flav, 0);
        for (const auto coupling : couplings) {
          if (coupling.first != flav) {
            auto ampl = std::make_pair(
                Cluster_Ampl_Key{{i, std::abs(coupling.first)}},
                CreateSU2TransformedAmplitude(
                    (*current_ampl), {{i, coupling.first}}));
            ampls.insert(std::move(ampl));
          }
        }
      }
    }

    if (activecoeffs.find(EWSudakov_Log_Type::lSSC) != activecoeffs.end()) {
      for (size_t k{0}; k < nlegs; ++k) {
        for (size_t l{0}; l < k; ++l) {
          // s-channel-related loops will have vanishing log coeffs
          if (k == 1 && l == 0)
            continue;
          if (nlegs == 4 && k == 3 && l == 2)
            continue;
          const auto kflav = (*current_ampl)->Leg(k)->Flav();
          const auto lflav = (*current_ampl)->Leg(l)->Flav();

          // I^Z * I^Z terms
          auto kcouplings = ewgroupconsts.IZ(kflav, 1);
          auto lcouplings = ewgroupconsts.IZ(lflav, 1);
            for (const auto kcoupling : kcouplings) {
              for (const auto lcoupling : lcouplings) {
                if (kcoupling.first != kflav ||
                    lcoupling.first != lflav) {
                  auto ampl = std::make_pair(
                      Cluster_Ampl_Key{{k, std::abs(kcoupling.first)},
                                       {l, std::abs(lcoupling.first)}},
                      CreateSU2TransformedAmplitude(
                          (*current_ampl), {{k, kcoupling.first}, {l, lcoupling.first}}));
                  ampls.insert(std::move(ampl));
                }
              }
            }

          // I^\pm * I^\pm terms
          for (size_t isplus{0}; isplus < 2; ++isplus) {
            kcouplings = ewgroupconsts.Ipm(kflav, 1, isplus);
            lcouplings = ewgroupconsts.Ipm(lflav, 1, !isplus);
            for (const auto kcoupling : kcouplings) {
              for (const auto lcoupling : lcouplings) {
                if (kcoupling.first != kflav ||
                    lcoupling.first != lflav) {
                  auto ampl = std::make_pair(
                      Cluster_Ampl_Key{{k, std::abs(kcoupling.first)},
                                       {l, std::abs(lcoupling.first)}},
                      CreateSU2TransformedAmplitude(
                          (*current_ampl), {{k, kcoupling.first}, {l, lcoupling.first}}));
                  ampls.insert(std::move(ampl));
                }
              }
            }
          }
        }
      }
    }
  }
  //for (auto& kv : ampls)
  //  Process_Base::SortFlavours(kv.second.get());
  return ampls;
}

Cluster_Amplitude_UP EWSudakov_Amplitudes::CreateAmplitude(Process_Base* proc)
{
  auto ampl = MakeClusterAmpl();
  ampl->SetNIn(proc->NIn());
  ampl->SetOrderQCD(proc->MaxOrder(0));
  for (size_t i(1);i<proc->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+proc->MaxOrder(i));
  for(int i(0);i<proc->NIn()+proc->NOut();++i)
    if (i<proc->NIn()) ampl->CreateLeg(Vec4D(),proc->Flavours()[i].Bar());
    else ampl->CreateLeg(Vec4D(),proc->Flavours()[i]);
  ampl->SetProc(proc);
  ampl->SetProcs(proc->AllProcs());
  return ampl;
}

Cluster_Amplitude_UP EWSudakov_Amplitudes::CreateSU2TransformedAmplitude(
    const Cluster_Amplitude_UP& ampl, const Leg_Kfcode_Map_Signed& flavs)
{
  auto campl = CopyClusterAmpl(ampl);
  for (const auto& kv : flavs) {
    campl->Leg(kv.first)->SetFlav(kv.second);
  }
  return campl;
}

EWSudakov_Amplitudes::Permutation_Map EWSudakov_Amplitudes::CreatePermutations(
    const Cluster_Amplitude_UPM& ampls)
{
  Permutation_Map permutations;
  for (const auto& kv : ampls) {
    std::vector<size_t> permutation;
    if (kv.first.empty()) {
      permutation.resize(kv.second->Legs().size());
      std::iota(std::begin(permutation), std::end(permutation), 0);
    } else {
      permutation = CalculateLegPermutation(kv.second);
    }
    permutations.insert(std::make_pair(kv.first, permutation));
  }
  return permutations;
}

std::vector<size_t> EWSudakov_Amplitudes::CalculateLegPermutation(
    const Cluster_Amplitude_UP& ampl)
{
  // TODO: This assumes that the base amplitude is ordered 0, 1, 2, ...
  std::vector<size_t> v;
  for (const auto* leg : ampl->Legs()) {
    if (IdCount(leg->Id()) > 1)
      // TODO: handle multi-ID legs
      THROW(not_implemented, "Clustering not supported yet.");
    v.push_back(ID(leg->Id()).front());
  }
  return v;
}
