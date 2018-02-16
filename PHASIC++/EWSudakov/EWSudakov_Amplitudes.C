#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

EWSudakov_Amplitudes::EWSudakov_Amplitudes(Process_Base* proc):
  ampls{ CreateAmplitudes(proc) },
  legpermutations{ CalculateLegPermutations(ampls) }
{
}

Cluster_Amplitude& EWSudakov_Amplitudes::Rotated(size_t legindex)
{
  const auto it = ampls.find(legindex);
  if (it == ampls.end())
    THROW(fatal_error, "Leg " + ToString(legindex) + "has not been rotated");
  return *ampls.find(legindex)->second;
}

std::vector<size_t>& EWSudakov_Amplitudes::LegPermutation(
    size_t legindex)
{
  const auto it = legpermutations.find(legindex);
  if (it == legpermutations.end())
    THROW(fatal_error, "Leg " + ToString(legindex) + "has not been rotated");
  return legpermutations.find(legindex)->second;
}

void EWSudakov_Amplitudes::UpdateMomenta(const ATOOLS::Vec4D_Vector& mom)
{
  for (auto& ampl : ampls) {
    for(int i(0); i < ampl.second->Legs().size(); ++i) {
      const auto j = (ampl.first == -1) ? i : LegPermutation(ampl.first)[i];
      ampl.second->Leg(i)->SetMom(mom[j]);
    }
  }
}

std::map<int, Cluster_Amplitude_UP> EWSudakov_Amplitudes::CreateAmplitudes(Process_Base* proc)
{
  std::map<int, Cluster_Amplitude_UP> ampls;
  auto itbool = ampls.insert(std::make_pair(-1, CreateAmplitude(proc)));
  auto& ampl = itbool.first->second;
  for (size_t i{ 0 }; i < ampl->Legs().size(); ++i) {
    const auto flav = ampl->Leg(i)->Flav();
    if (flav.IsPhoton() || flav.Kfcode() == kf_Z)
      // TODO: include more flavours when needed
      ampls.insert(std::make_pair(i, CreateSU2RotatedAmplitude(ampl, i)));
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

Cluster_Amplitude_UP EWSudakov_Amplitudes::CreateSU2RotatedAmplitude(
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

std::map<size_t, std::vector<size_t>> EWSudakov_Amplitudes::CalculateLegPermutations(
    const std::map<int, Cluster_Amplitude_UP>& ampls)
{
  std::map<size_t, std::vector<size_t>> permutations;
  for (const auto& kv : ampls) {
    if (kv.first == -1)
      continue;
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
