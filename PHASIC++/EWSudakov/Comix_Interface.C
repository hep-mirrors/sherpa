#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;
using namespace MODEL;

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

void Comix_Interface::FillSpinAmplitudes(
    std::vector<Spin_Amplitudes>& spinampls,
    const ATOOLS::Cluster_Amplitude& ampl) const
{
  Cluster_Amplitude_UP campl {CopyClusterAmpl(ampl)};
  // TODO: can't we set these scales to the values used for the original
  // calculation (pre-Kfactor)? Is there any reason we set them to Ecms? Does
  // it matter?
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  Process_Base* proc = GetProcess(*campl);
  if (proc == nullptr)
    return;
  proc->Differential(*campl, Weight_Type::nominal, m_differentialmode);
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
    if (Flavour(kf_jet).Includes(fl))
      fl = Flavour(kf_jet);
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl, "", ""));
  }
  for (size_t i{ampl->NIn()}; i < ampl->Legs().size(); ++i) {
    Flavour fl{ampl->Leg(i)->Flav()};
    if (Flavour(kf_jet).Includes(fl))
      fl = Flavour(kf_jet);
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl, "", ""));
  }

  // copy coupling orders and allow for any SMGold coupling order
  pi.m_maxcpl = p_proc->Info().m_maxcpl;
  pi.m_mincpl = p_proc->Info().m_mincpl;
  pi.m_maxcpl.push_back(99);
  pi.m_mincpl.push_back(0);

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
