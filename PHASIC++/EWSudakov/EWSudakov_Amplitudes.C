#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

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
EWSudakov_Amplitudes::BaseAmplitude(std::vector<int> spincombination) noexcept
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
  // TODO: also, reuse existing implementation elsewhere if existing
  assert(BaseAmplitude().NIn() == 2 && BaseAmplitude().Legs().size() == 4);
  const auto Etot2 = MandelstamS();
  for (auto& ampl : ampls) {
    if (ampl.first == s_baseamplkey)
      continue;
    const auto& permutation = LegPermutation(ampl.first);

    // set incoming momenta
    // TODO: drop the implicit assumption that all initial-state particles are
    // massless (?)
    for (int i{ 0 }; i < ampl.second->NIn(); ++i) {
      ampl.second->Leg(i)->SetMom(-mom[permutation[i]]);
    }

    // calc outgoing momenta such that all particles are on their mass shell;
    // to guarantee this we modify the energy and the absolute value of the
    // momenta, but not the directions of the momenta
    const Vec3D out_mom {BaseAmplitude().Leg(permutation[2])->Mom()};
    const Vec3D normed_out_mom {out_mom / out_mom.Abs()};
    const auto m22 = sqr(ATOOLS::Abs<double>(ampl.second->Leg(2)->Flav().Mass()));
    const auto m32 = sqr(ATOOLS::Abs<double>(ampl.second->Leg(3)->Flav().Mass()));
    // calc the fraction of energy that will be assigned to leg 2
    const auto x = (Etot2 + m22 - m32) / (2*Etot2);
    // calc the absolute momentum assigned to leg 2
    const auto p
      = sqrt((Etot2*Etot2 - 2*Etot2*(m22 + m32) + sqr(m22 - m32)) / (4*Etot2));
    ampl.second->Leg(2)->SetMom(Vec4D{x*sqrt(Etot2), p*normed_out_mom});
    ampl.second->Leg(3)->SetMom(Vec4D{(1-x)*sqrt(Etot2), -p*normed_out_mom});
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

  // create unmodified amplitude
  const auto& baseampl
    = ampls.insert(std::make_pair(s_baseamplkey,
                                  CreateAmplitude(proc))).first->second;
  const auto nlegs = baseampl->Legs().size();

  // create amplitudes where W and Z are replaced with their Goldstone
  // partners, phi and chi

  // store W and Z indizes
  std::vector<size_t> bosonindexes;
  bosonindexes.reserve(nlegs);
  for (size_t k{0}; k < nlegs; ++k) {
    const auto flav = baseampl->Leg(k)->Flav();
    if (flav.Kfcode() == kf_Z || flav.Kfcode() == kf_Wplus) {
      bosonindexes.push_back(k);
    }
  }
  // permute over replacing / not replacing each boson index, leaving out the
  // "do not replace anything" permutation (i.e. do not start with 1)
  size_t first_invalid_permutation{
      static_cast<size_t>(1 << bosonindexes.size())};
  for (size_t permutation{1}; permutation != first_invalid_permutation;
       ++permutation) {
    Leg_Kfcode_Map leg_set;
    Leg_Kfcode_Map_Signed leg_set_signed;
    for (size_t k{0}; k < nlegs; ++k) {
      if (permutation & (1 << k)) {
        const auto& flav = baseampl->Leg(k)->Flav();
        const auto flavcode = (long int)flav.GoldstoneBosonPartner();
        leg_set.emplace(k, std::abs(flavcode));
        leg_set_signed.emplace(k, flavcode);
      }
    }
    const auto it = ampls.find(leg_set);
    if (it == ampls.end()) {
      auto ampl = std::make_pair(
          leg_set, CreateSU2TransformedAmplitude(baseampl, leg_set_signed));
      ampls.insert(std::move(ampl));
    }
  }

  // create amplitudes needed for Z/photon interference terms
  if (activecoeffs.find(EWSudakov_Log_Type::Ls) != activecoeffs.end()) {
    for (size_t i{ 0 }; i < nlegs; ++i) {
      const auto flav = baseampl->Leg(i)->Flav();
      int newkf{ kf_none };
      if (flav.IsPhoton())
        newkf = kf_Z;
      else if (flav.Kfcode() == kf_Z)
        newkf = kf_photon;
      if (newkf != kf_none) {
        auto ampl =
            std::make_pair(Cluster_Ampl_Key{{i, newkf}},
                           CreateSU2TransformedAmplitude(
                               baseampl, {{i, static_cast<long int>(newkf)}}));
        ampls.insert(std::move(ampl));
      }
    }
  }

  // create amplitudes needed for W loops
  if (activecoeffs.find(EWSudakov_Log_Type::lSSC) != activecoeffs.end()) {
    for (size_t k{ 0 }; k < nlegs; ++k) {
      for (size_t l{ 0 }; l < k; ++l) {

        // s-channel-related loops will have vanishing log coeffs
        if (k == 1 && l == 0)
          continue;
        if (nlegs == 4 && k == 3 && l == 2)
          continue;

        const auto kflav = baseampl->Leg(k)->Flav();
        const auto lflav = baseampl->Leg(l)->Flav();

        // collect transformed particles
        Flavour_Vector newkflavs, newlflavs;
        if (kflav.IsLepton() || kflav.IsQuark()) {
          newkflavs.push_back(kflav.IsoWeakPartner());
        } else if (kflav.Kfcode() == kf_Wplus) {
          newkflavs.push_back(Flavour(kf_Z));  // for any W
          newkflavs.push_back(Flavour(kf_photon));  // for transverse W
          newkflavs.push_back(Flavour(kf_h0));  // for longitudinal W
        } else if (kflav.Kfcode() == kf_Z || kflav.Kfcode() == kf_photon) {
          newkflavs.push_back(Flavour(kf_Wplus, false));
          newkflavs.push_back(Flavour(kf_Wplus, true));
        }
        if (lflav.IsLepton() || lflav.IsQuark()) {
          newlflavs.push_back(lflav.IsoWeakPartner());
        } else if (lflav.Kfcode() == kf_Wplus) {
          newlflavs.push_back(Flavour(kf_Z));  // for any W
          newlflavs.push_back(Flavour(kf_photon));  // for transverse W
          newlflavs.push_back(Flavour(kf_h0));  // for longitudinal W
        } else if (kflav.Kfcode() == kf_Z || lflav.Kfcode() == kf_photon) {
          newlflavs.push_back(Flavour(kf_Wplus, false));
          newlflavs.push_back(Flavour(kf_Wplus, true));
        }

        // create valid amplitudes
        for (const auto newkflav : newkflavs) {
          for (const auto newlflav : newlflavs) {
            if (kflav.IntCharge() + lflav.IntCharge()
                != newkflav.IntCharge() + newlflav.IntCharge())
              continue;
            auto ampl =
                std::make_pair(Cluster_Ampl_Key{{k, newkflav.Kfcode()},
                                                {l, newlflav.Kfcode()}},
                               CreateSU2TransformedAmplitude(
                                   baseampl, {{k, newkflav}, {l, newlflav}}));
            ampls.insert(std::move(ampl));
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
