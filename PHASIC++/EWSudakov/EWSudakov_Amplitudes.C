#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/EWSudakov/EWGroupConstants.H"
#include "PHASIC++/Process/Process_Base.H"

#include <numeric>
#include <cassert>

using namespace PHASIC;
using namespace ATOOLS;

const Cluster_Ampl_Key EWSudakov_Amplitudes::s_baseamplkey
= Leg_Kfcode_Map{};

EWSudakov_Amplitudes::EWSudakov_Amplitudes(
    Process_Base* proc, const std::set<EWSudakov_Log_Type>& activecoeffs)
    : ampls {CreateAmplitudes(proc, activecoeffs)}
{
  // fill helper maps that are used to expose specific sets of amplitudes
  for (const auto& kv : ampls) {
    all_ampls[kv.first] = kv.second.get();
    bool is_goldstone_only {true};
    for (const auto& leg_kfcode_pair : kv.first) {
      if (!(leg_kfcode_pair.second == kf_phiplus ||
          leg_kfcode_pair.second == kf_chi)) {
        is_goldstone_only = false;
        break;
      }
    }
    if (is_goldstone_only) {
      goldstone_only_ampls[kv.first] = kv.second.get();
    }
  }
}

EWSudakov_Amplitudes::~EWSudakov_Amplitudes()
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::BaseAmplitude() noexcept
{
  return SU2TransformedAmplitude(s_baseamplkey);
}

Cluster_Amplitude&
EWSudakov_Amplitudes::BaseAmplitude(std::vector<int> spincombination)
{
  return SU2TransformedAmplitude(GoldstoneBosonReplacements(spincombination));
}

Leg_Kfcode_Map EWSudakov_Amplitudes::GoldstoneBosonReplacements(
    std::vector<int> spincombination)
{
  Leg_Kfcode_Map leg_set;
  for (int i {0}; i < NumberOfLegs(); ++i) {
    if (spincombination[i] == 2) {
      const auto flav = BaseAmplitude().Leg(i)->Flav();
      if (flav.Kfcode() == kf_Z || flav.Kfcode() == kf_Wplus) {
        leg_set.emplace(i, flav.GoldstoneBosonPartner().Kfcode());
      }
    }
  }
  return leg_set;
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

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  DEBUG_FUNC(mom);
  for (int i {0}; i < NumberOfLegs(); ++i) {
    BaseAmplitude().SetMom(i, mom[i]);
  }

  particles.reserve(NumberOfLegs());
  for (int n {0}; n < NumberOfLegs(); n++)
    particles.push_back(new Particle{});

  for (auto& ampl : ampls) {
    if (ampl.first == s_baseamplkey)
      continue;

    Vec4D_Vector new_moms;
    new_moms.reserve(NumberOfLegs());
    for (int j {0}; j < NumberOfLegs(); ++j) {
      Vec4D mom = BaseAmplitude().Mom(j);
      particles[j]->SetFinalMass(ampl.second->Leg(j)->Flav().Mass());
      new_moms.push_back(mom);
    }
    const auto did_stretch = stretcher.StretchMomenta(particles, new_moms);
    if (!did_stretch) {
      ampl.second->SetFlag(1 << 4);
      continue;
    }
    for (int j {0}; j < NumberOfLegs(); ++j) {
      ampl.second->SetMom(j, new_moms[j]);
    }
  }

  for (auto* p : particles)
    delete p;
  particles.clear();
}

void EWSudakov_Amplitudes::UpdateColors(const Int_Vector& I, const Int_Vector& J)
{
  for (auto& ampl : ampls) {
    for (int j {0}; j < NumberOfLegs(); ++j) {
      ampl.second->Leg(j)->SetCol(ColorID(I[j], J[j]));
    }
  }
}

double EWSudakov_Amplitudes::MandelstamS()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Mom(0) + ampl.Mom(1)).Abs2();
}

double EWSudakov_Amplitudes::MandelstamT()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Mom(0) - ampl.Mom(2)).Abs2();
}

double EWSudakov_Amplitudes::MandelstamU()
{
  const auto& ampl = BaseAmplitude();
  return (ampl.Mom(0) - ampl.Mom(3)).Abs2();
}

Cluster_Amplitude_UPM EWSudakov_Amplitudes::CreateAmplitudes(
    Process_Base* proc, const std::set<EWSudakov_Log_Type>& activecoeffs) const
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
    auto* current_ampl = &baseampl;
    Leg_Kfcode_Map goldstone_leg_set;
    if (permutation != 0) {
      Leg_Kfcode_Map_Signed goldstone_leg_set_signed;
      for (size_t k{0}; k < bosonindexes.size(); ++k) {
        if (permutation & (1 << k)) {
          const auto bosonindex = bosonindexes[k];
          const auto& flav = (*current_ampl)->Leg(bosonindex)->Flav();
          const auto flavcode = (long int)flav.GoldstoneBosonPartner();
          goldstone_leg_set.emplace(bosonindex, std::abs(flavcode));
          goldstone_leg_set_signed.emplace(bosonindex, flavcode);
        }
      }
      auto it = ampls.find(goldstone_leg_set);
      if (it == ampls.end()) {
        auto ampl = std::make_pair(
            goldstone_leg_set, CreateSU2TransformedAmplitude((*current_ampl), goldstone_leg_set_signed));
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
          auto leg_set = Cluster_Ampl_Key {{i, newkf}};
          leg_set.insert(std::begin(goldstone_leg_set),
                         std::end(goldstone_leg_set));
          auto ampl = std::make_pair(
              leg_set,
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
            auto leg_set = Cluster_Ampl_Key {{i, std::abs(coupling.first)}};
            leg_set.insert(std::begin(goldstone_leg_set),
                           std::end(goldstone_leg_set));
            auto ampl = std::make_pair(
                leg_set,
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
                  auto leg_set =
                      Cluster_Ampl_Key {{k, std::abs(kcoupling.first)},
                                        {l, std::abs(lcoupling.first)}};
                  leg_set.insert(std::begin(goldstone_leg_set),
                                 std::end(goldstone_leg_set));
                  auto ampl = std::make_pair(
                      leg_set,
                      CreateSU2TransformedAmplitude(
                          (*current_ampl),
                          {{k, kcoupling.first}, {l, lcoupling.first}}));
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
                  auto leg_set =
                      Cluster_Ampl_Key {{k, std::abs(kcoupling.first)},
                                        {l, std::abs(lcoupling.first)}};
                  leg_set.insert(std::begin(goldstone_leg_set),
                                 std::end(goldstone_leg_set));
                  auto ampl = std::make_pair(
                      leg_set,
                      CreateSU2TransformedAmplitude(
                          (*current_ampl),
                          {{k, kcoupling.first}, {l, lcoupling.first}}));
                  ampls.insert(std::move(ampl));
                }
              }
            }
          }
        }
      }
    }
  }
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
