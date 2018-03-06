#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <numeric>

using namespace PHASIC;
using namespace ATOOLS;

const EWSudakov_Amplitudes::Cluster_Ampl_Key EWSudakov_Amplitudes::s_baseamplkey
= std::make_pair(EWSudakov_Amplitude_Type::Base, Leg_Index_Set{});

EWSudakov_Amplitudes::EWSudakov_Amplitudes(Process_Base* proc):
  ampls{ CreateAmplitudes(proc) },
  permutations{ CreatePermutations(ampls) }
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::Unrotated() noexcept
{
  return Rotated(s_baseamplkey.first, s_baseamplkey.second);
}

Cluster_Amplitude& EWSudakov_Amplitudes::Rotated(EWSudakov_Amplitude_Type type,
                                                 Leg_Index_Set indizes)
{
  const auto key = std::make_pair(type, indizes);
  const auto it = ampls.find(key);
  if (it == ampls.end())
    THROW(fatal_error, "Rotated amplitude not found");
  return *(it->second);
}

std::vector<size_t>& EWSudakov_Amplitudes::LegPermutation(
    EWSudakov_Amplitude_Type type, Leg_Index_Set indizes)
{
  const auto key = std::make_pair(type, indizes);
  const auto it = permutations.find(key);
  if (it == permutations.end())
    THROW(fatal_error, "Permutation not found");
  return it->second;
}

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  for (auto& ampl : ampls) {
    const auto ampltype = ampl.first.first;
    const auto legindizes = ampl.first.second;
    const auto& permutation = LegPermutation(ampltype, legindizes);
    for(int i{ 0 }; i < ampl.second->Legs().size(); ++i) {
      ampl.second->Leg(i)->SetMom(mom[permutation[i]]);
    }
  }
}

EWSudakov_Amplitudes::Cluster_Amplitude_UPM
EWSudakov_Amplitudes::CreateAmplitudes(Process_Base* proc) const
{
  Cluster_Amplitude_UPM ampls;
  const auto& baseampl
    = ampls.insert(std::make_pair(s_baseamplkey, CreateAmplitude(proc))).first->second;
  // create amplitudes needed for Z/photon interference terms
  for (size_t i{ 0 }; i < baseampl->Legs().size(); ++i) {
    const auto flav = baseampl->Leg(i)->Flav();
    if (flav.IsPhoton() || flav.Kfcode() == kf_Z) {
      const auto key
        = std::make_pair(EWSudakov_Amplitude_Type::LSCZ, Leg_Index_Set{i});
      ampls.insert(std::make_pair(key, CreateLSCZAmplitude(baseampl, i)));
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

Cluster_Amplitude_UP EWSudakov_Amplitudes::CreateLSCZAmplitude(
    const Cluster_Amplitude_UP& ampl, size_t legindex)
{
  auto campl = CopyClusterAmpl(ampl);
  auto* leg = campl->Leg(legindex);
  auto flav = leg->Flav();
  Flavour newflav;
  // TODO: generalise to other flavours, and, when generalised, migrate to
  // Flavour class where it belongs
  if (flav.IsPhoton())
    newflav = Flavour{kf_Z};
  else if (flav.Kfcode() == kf_Z)
    newflav = Flavour{kf_photon};
  leg->SetFlav(newflav);
  return campl;
}

EWSudakov_Amplitudes::Permutation_Map EWSudakov_Amplitudes::CreatePermutations(
    const Cluster_Amplitude_UPM& ampls)
{
  Permutation_Map permutations;
  for (const auto& kv : ampls) {
    std::vector<size_t> permutation;
    switch (kv.first.first) {
      case EWSudakov_Amplitude_Type::Base:
        permutation.resize(kv.second->Legs().size());
        std::iota(std::begin(permutation), std::end(permutation), 0);
        break;
      case EWSudakov_Amplitude_Type::LSCZ:
        permutation = CalculateLSCZLegPermutation(kv.second);
        break;
      default:
        THROW(not_implemented, "Missing implementation");
    }
    permutations.insert(std::make_pair(kv.first, permutation));
  }
  return permutations;
}

std::vector<size_t> EWSudakov_Amplitudes::CalculateLSCZLegPermutation(
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
