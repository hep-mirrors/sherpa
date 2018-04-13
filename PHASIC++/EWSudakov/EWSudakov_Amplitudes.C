#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <numeric>
#include <cassert>

using namespace PHASIC;
using namespace ATOOLS;

const EWSudakov_Amplitudes::Cluster_Ampl_Key EWSudakov_Amplitudes::s_baseamplkey
= Leg_Set{};

EWSudakov_Amplitudes::EWSudakov_Amplitudes(
    Process_Base* proc, const std::set<EWSudakov_Log_Type>& activecoeffs):
  ampls{ CreateAmplitudes(proc, activecoeffs) },
  permutations{ CreatePermutations(ampls) }
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::Unrotated() noexcept
{
  return Rotated(s_baseamplkey);
}

Cluster_Amplitude& EWSudakov_Amplitudes::Rotated(const Leg_Set& legs)
{
  const auto it = ampls.find(legs);
  if (it == ampls.end())
    THROW(fatal_error, "Rotated amplitude not found");
  return *(it->second);
}

std::vector<size_t>& EWSudakov_Amplitudes::LegPermutation(const Leg_Set& legs)
{
  const auto it = permutations.find(legs);
  if (it == permutations.end())
    THROW(fatal_error, "Permutation not found");
  return it->second;
}

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  for (auto& ampl : ampls) {
    const auto& permutation = LegPermutation(ampl.first);
    for(int i{ 0 }; i < ampl.second->Legs().size(); ++i) {
      ampl.second->Leg(i)->SetMom(mom[permutation[i]]);
    }
  }
}

double EWSudakov_Amplitudes::MandelstamS()
{
  const auto& ampl = Unrotated();
  return (ampl.Leg(0)->Mom() + ampl.Leg(1)->Mom()).Abs2();
}

double EWSudakov_Amplitudes::MandelstamT()
{
  const auto& ampl = Unrotated();
  return (ampl.Leg(0)->Mom() - ampl.Leg(2)->Mom()).Abs2();
}

double EWSudakov_Amplitudes::MandelstamU()
{
  const auto& ampl = Unrotated();
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

  // create amplitudes needed for Z/photon interference terms
  if (activecoeffs.find(EWSudakov_Log_Type::Ls) != activecoeffs.end()) {
    for (size_t i{ 0 }; i < baseampl->Legs().size(); ++i) {
      const auto flav = baseampl->Leg(i)->Flav();
      int newkf{ kf_none };
      if (flav.IsPhoton())
        newkf = kf_Z;
      else if (flav.Kfcode() == kf_Z)
        newkf = kf_photon;
      if (newkf != kf_none) {
        auto ampl = std::make_pair(
            Leg_Set{ {i, newkf} },
            CreateRotatedAmplitude(baseampl, {{i, Flavour(newkf)}}));
        ampls.insert(std::move(ampl));
      }
    }
  }

  // create amplitudes needed for W loops
  if (activecoeffs.find(EWSudakov_Log_Type::lSSC) != activecoeffs.end()) {
    for (size_t k{ 0 }; k < baseampl->Legs().size(); ++k) {
      for (size_t l{ 0 }; l < k; ++l) {
        const auto kflav = baseampl->Leg(k)->Flav();
        const auto lflav = baseampl->Leg(l)->Flav();

        // collect transformed particles
        Flavour_Vector newkflavs, newlflavs;
        if (kflav.IsLepton() || kflav.IsQuark()) {
          newkflavs.push_back(kflav.IsoWeakPartner());
        } else if (kflav.Kfcode() == kf_Wplus) {
          newkflavs.push_back(Flavour(kf_Z));
          newkflavs.push_back(Flavour(kf_photon));
        }
        if (lflav.IsLepton() || lflav.IsQuark()) {
          newlflavs.push_back(lflav.IsoWeakPartner());
        } else if (lflav.Kfcode() == kf_Wplus) {
          newlflavs.push_back(Flavour(kf_Z));
          newlflavs.push_back(Flavour(kf_photon));
        }

        // create valid amplitudes
        for (const auto newkflav : newkflavs) {
          for (const auto newlflav : newlflavs) {
            if (kflav.IntCharge() + lflav.IntCharge()
                != newkflav.IntCharge() + newlflav.IntCharge())
              continue;
	    auto ampl = std::make_pair(
                Leg_Set{ {k, newkflav.Kfcode()}, {l, newlflav.Kfcode()} },
                CreateRotatedAmplitude(baseampl, {{k, newkflav}, {l, newlflav}}));
            ampls.insert(std::move(ampl));
          }
        }
      }
    }
  }
  for (auto& kv : ampls)
    Process_Base::SortFlavours(kv.second.get());
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

Cluster_Amplitude_UP
EWSudakov_Amplitudes::CreateRotatedAmplitude(const Cluster_Amplitude_UP& ampl,
                                             const LegIndex_Flavour_Map& flavs)
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
  // TODO: This assumes that the unrotated amplitude is ordered 0, 1, 2, ...
  std::vector<size_t> v;
  for (const auto* leg : ampl->Legs()) {
    if (IdCount(leg->Id()) > 1)
      // TODO: handle multi-ID legs
      THROW(not_implemented, "Clustering not supported yet.");
    v.push_back(ID(leg->Id()).front());
  }
  return v;
}
