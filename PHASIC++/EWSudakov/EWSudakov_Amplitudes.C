#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

EWSudakov_Amplitudes::EWSudakov_Amplitudes(Process_Base* proc):
  baseampl{ CreateAmplitude(proc) },
  zphotoninterferenceampls{ CreateZPhotonInterferenceAmplitudes(proc) },
  zphotonlegpermutations{ CalculateLegPermutations(zphotoninterferenceampls) }
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::Rotated(EWSudakov_Amplitude_Type type,
                                                 size_t legindex)
{
  switch (type) {
    case EWSudakov_Amplitude_Type::Base:
      return Unrotated();
    case EWSudakov_Amplitude_Type::ZPhotonInterference: {
      const auto it = zphotoninterferenceampls.find(legindex);
      if (it == zphotoninterferenceampls.end())
        THROW(fatal_error, "Rotated amplitude not found");
      return *(it->second);
    }
    default:
      THROW(not_implemented, "EWSudakov Amplitude type not implemented");
  }
}

std::vector<size_t>& EWSudakov_Amplitudes::LegPermutation(
    EWSudakov_Amplitude_Type type, size_t legindex)
{
  switch (type) {
    case EWSudakov_Amplitude_Type::ZPhotonInterference: {
      const auto it = zphotonlegpermutations.find(legindex);
      if (it == zphotonlegpermutations.end())
        THROW(fatal_error, "Permutation not found");
      return it->second;
    }
    default:
      THROW(not_implemented, "EWSudakov Amplitude type not implemented");
  }
}

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  for(int i{ 0 }; i < baseampl->Legs().size(); ++i) {
    baseampl->Leg(i)->SetMom(mom[i]);
  }
  for (auto& ampl : zphotoninterferenceampls) {
    const auto& permutation =
      LegPermutation(EWSudakov_Amplitude_Type::ZPhotonInterference, ampl.first);
    for(int i{ 0 }; i < ampl.second->Legs().size(); ++i) {
      ampl.second->Leg(i)->SetMom(mom[permutation[i]]);
    }
  }
}

EWSudakov_Amplitudes::Cluster_Amplitude_UPM
EWSudakov_Amplitudes::CreateZPhotonInterferenceAmplitudes(Process_Base* proc) const
{
  Cluster_Amplitude_UPM ampls;
  // create amplitudes needed for Z/photon interference terms
  for (size_t i{ 0 }; i < baseampl->Legs().size(); ++i) {
    const auto flav = baseampl->Leg(i)->Flav();
    if (flav.IsPhoton() || flav.Kfcode() == kf_Z) {
      ampls.insert(
          std::make_pair(i, CreateZPhotonInterferenceAmplitude(baseampl, i)));
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

Cluster_Amplitude_UP EWSudakov_Amplitudes::CreateZPhotonInterferenceAmplitude(
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

std::map<size_t, std::vector<size_t>>
EWSudakov_Amplitudes::CalculateLegPermutations(
    const Cluster_Amplitude_UPM& ampls)
{
  std::map<size_t, std::vector<size_t>> permutations;
  for (const auto& kv : ampls) {
    permutations[kv.first] = CalculateLegPermutation(kv.second);
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

ClusterAmplitude_Vector EWSudakov_Amplitudes::AllAmplitudes() noexcept
{
  ClusterAmplitude_Vector v;
  v.push_back(baseampl.get());
  for (const auto& ampl : zphotoninterferenceampls)
    v.push_back(ampl.second.get());
  return v;
}
