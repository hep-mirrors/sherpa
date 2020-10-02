#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "COMIX/Main/Single_Process.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;
using namespace MODEL;
using namespace METOOLS;

NLOTypeStringProcessMap_Map Comix_Interface::s_apmap;

Comix_Interface::Comix_Interface(Process_Base* proc,
                                 const EWSudakov_Amplitudes& ampls)
    : p_proc {proc}, m_procname_suffix {"Sudakov"}, m_differentialmode {2 | 4}
{
  AdaptToProcessColorScheme();
  InitializeProcesses(ampls.All());
}

Comix_Interface::Comix_Interface(Process_Base* proc,
                                 const std::string& procname_suffix)
    : p_proc {proc},
      m_procname_suffix {procname_suffix},
      m_differentialmode {2 | 4}
{
  AdaptToProcessColorScheme();
}

Complex Comix_Interface::GetSpinAmplitude(ATOOLS::Cluster_Amplitude& ampl,
                                          const std::vector<int>& spins) const
{
  Cluster_Amplitude_UP campl {CopyClusterAmpl(ampl)};

  // during process set-up, flavours (in the proc info) are implicitly sorted,
  // hence we also need to sort them here, where we step out of the EWSudakov
  // module (in which we don't order flavours, to have a simpler book-keeping
  Process_Base::SortFlavours(campl.get());

  auto ait = m_spinampls.find(&ampl);
  if (ait == m_spinampls.end()) {
    FillSpinAmplitudes(m_spinampls[&ampl], *campl);
    ait = m_spinampls.find(&ampl);
  }
  if (ait->second.empty())
    return 0.0;
  assert(ait->second.size() == 1);

  auto pit = m_permutations.find(&ampl);
  if (pit == m_permutations.end()) {
    // calculate leg permutation (in particular by calling SortFlavours above)
    for (const auto* leg : campl->Legs()) {
      if (IdCount(leg->Id()) > 1)
        // TODO: handle multi-ID legs
        THROW(not_implemented, "Clustering not supported yet.");
      m_permutations[&ampl].push_back(ID(leg->Id()).front());
    }
    pit = m_permutations.find(&ampl);
  }

  // apply the permutation to the spin combination
  std::vector<int> permutatedspins;
  for (const auto& idx : pit->second) {
    permutatedspins.push_back(spins[idx]);
  }

  auto sit = m_permutation_signs.find(&ampl);
  if (sit == m_permutation_signs.end()) {
    // calculate the sign of the permutation
    auto sign = Comix_Interface::CalcPermutationSign(pit->second, ampl);
    m_permutation_signs[&ampl] = sign;
    sit = m_permutation_signs.find(&ampl);
  }

  // retrieve correct ME entry
  return sit->second * ait->second[0].Get(permutatedspins);
}

void Comix_Interface::FillSpinAmplitudes(std::vector<Spin_Amplitudes>& spinampls,
                                         const ATOOLS::Cluster_Amplitude& ampl) const
{
  Cluster_Amplitude_UP campl {CopyClusterAmpl(ampl)};
  SetScales(*campl);
  Process_Base* proc = GetProcess(*campl);
  if (proc == nullptr)
    return;
  proc->Differential(*campl, Variations_Mode::nominal_only, m_differentialmode);
  std::vector<std::vector<Complex>> cols;
  proc->FillAmplitudes(spinampls, cols);
}

PHASIC::Process_Base*
Comix_Interface::GetProcess(const ATOOLS::Cluster_Amplitude& ampl) const
{
  const auto loprocmapit = ProcessMap().find(nlo_type::lo);
  if (loprocmapit == ProcessMap().end())
    return nullptr;
  std::string pname {Process_Base::GenerateName(&ampl)};
  auto pit = loprocmapit->second->find(pname);
  if (pit == loprocmapit->second->end()) {
    return nullptr;
  }
  return pit->second;
}

void Comix_Interface::InitializeProcesses(const Cluster_Amplitude_PM& ampls)
{
  DEBUG_FUNC("");
  auto& s = Settings::GetMainSettings();
  const auto graph_path =
      s["PRINT_EWSUDAKOV_GRAPHS"].SetDefault("").Get<std::string>();
  for (const auto& kv : ampls) {
    const auto& ampl = kv.second;
    PHASIC::Process_Base* proc = GetProcess(*ampl);
    if (proc != nullptr)
      // ignore processes that have already been initialized
      continue;
    const Process_Info pi = CreateProcessInfo(ampl, graph_path);
    InitializeProcess(pi);
  }
}

Process_Info
Comix_Interface::CreateProcessInfo(const Cluster_Amplitude* ampl,
                                   const std ::string& graph_path)
{
  Process_Info pi;
  pi.m_addname = "__" + m_procname_suffix;
  pi.m_megenerator = "Comix";
  if (graph_path != "") {
    pi.m_gpath = graph_path;
  }

  // set external particles
  // TODO: do we need the parton -> jet transformation?
  for (size_t i{0}; i < ampl->NIn(); ++i) {
    Flavour fl(ampl->Leg(i)->Flav().Bar());
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl, "", ""));
  }
  for (size_t i{ampl->NIn()}; i < ampl->Legs().size(); ++i) {
    Flavour fl{ampl->Leg(i)->Flav()};
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl, "", ""));
  }

  // copy coupling orders and allow for any SMGold coupling order
  pi.m_maxcpl = p_proc->Info().m_maxcpl;
  pi.m_mincpl = p_proc->Info().m_mincpl;
  pi.m_maxacpl = p_proc->Info().m_maxacpl;
  pi.m_minacpl = p_proc->Info().m_minacpl;
  pi.m_maxcpl[2] = 99;
  pi.m_mincpl[2] = 0;
  pi.m_maxacpl[2] = 99;
  pi.m_minacpl[2] = 0;

  return pi;
}

void Comix_Interface::InitializeProcess(const Process_Info& pi)
{
  auto proc = p_proc->Generator()->Generators()->InitializeProcess(pi, false);
  if (proc == NULL) {
    msg_Error() << "WARNING: Comix_Interface::InitializeProcess can not"
                << "initialize process for process info: " << pi << '\n';
    return;
  }
  proc->SetSelector(Selector_Key{});
  proc->SetScale(Scale_Setter_Arguments(
      MODEL::s_model, "VAR{" + ToString(sqr(rpa->gen.Ecms())) + "}",
      "Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  //proc->Get<COMIX::Process_Base>()->Tests();
  proc->FillProcessMap(&ProcessMap());
}

void Comix_Interface::AdaptToProcessColorScheme()
{
  if (p_proc->Integrator()->ColorScheme() == cls::sum) {
    // for some reason, this is needed when summing colours; if colour are
    // sampled, however, we can't have this because it then triggers the
    // generation of a new random colour point each time we call Differential()
    m_differentialmode |= 128;
  }
}

void Comix_Interface::SetScales(ATOOLS::Cluster_Amplitude& ampl) const
{
  // TODO: can't we set these scales to the values used for the original
  // calculation (pre-Kfactor)? Is there any reason we set them to Ecms? Does
  // it matter?
  const auto scale2 = sqr(rpa->gen.Ecms());
  ampl.SetMuR2(scale2);
  ampl.SetMuF2(scale2);
  ampl.SetMuQ2(scale2);
}

double Comix_Interface::CalcPermutationSign(std::vector<size_t> perm,
    const Cluster_Amplitude& ampl)
{
  // extract permutation of final-state fermions only
  static const size_t boson_idx = std::numeric_limits<size_t>::max();
  perm.erase(perm.begin(), perm.begin()+ampl.NIn());
  for (auto i = ampl.NIn(); i < ampl.Legs().size(); ++i) {
    if (ampl.Leg(i)->Flav().IsBoson()) {
      auto idx = ID(ampl.Leg(i)->Id()).front();
      for (auto it = perm.begin(); it != perm.end(); ++it) {
        if (*it == idx) {
          *it = boson_idx;
        } else if (*it > idx && *it != boson_idx) {
          --(*it);
        }
      }
    }
  }
  for (auto& p : perm) {
    if (p != boson_idx)
      p -= ampl.NIn();
  }
  for (auto it = perm.begin(); it != perm.end(); ) {
    if (*it == boson_idx)
      it = perm.erase(it);
    else
      ++it;
  }

  // calculate number of even-length cycles
  std::vector<bool> checklist(perm.size(), false);
  size_t number_of_even_cycles {0};
  for (int i {0}; i < perm.size(); ++i) {
    if (checklist[i])
      continue;
    checklist[i] = true;
    auto j = i;
    int cycle_length = 1;
    while (!checklist[perm[j]]) {
      j = perm[j];
      cycle_length++;
    }
    if (cycle_length % 2 == 0)
      ++number_of_even_cycles;
  }

  // return sign of the fermion permutation
  if (number_of_even_cycles % 2 == 1)
    return -1.0;
  else
    return 1.0;
}
