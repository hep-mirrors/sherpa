#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "ATOOLS/Phys/Color.H"
#include "COMIX/Main/Single_Process.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Shell_Tools.H"

#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;
using namespace MODEL;

Comix_Interface::Comix_Interface(Process_Base* proc,
                                 const EWSudakov_Amplitudes& ampls)
    : p_proc{proc}, m_procname_suffix{"Sudakov"}
{
  InitializeProcesses(ampls.All());
}

Comix_Interface::Comix_Interface(Process_Base* proc,
                                 const std::string& procname_suffix)
    : p_proc {proc}, m_procname_suffix {procname_suffix}
{
}

void Comix_Interface::FillSpinAmplitudes(
    std::vector<Spin_Amplitudes>& spinampls,
    ATOOLS::Cluster_Amplitude& ampl) const
{
  const auto loprocmapit = m_apmap.find(nlo_type::lo);
  if (loprocmapit == m_apmap.end())
    THROW(fatal_error, "LO entry in process map not found");
  Cluster_Amplitude* campl(ampl.Copy());
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  std::string pname(Process_Base::GenerateName(campl));
  auto pit = loprocmapit->second->find(pname);
  if (pit->second == NULL) {
    msg_Error() << "Looking for amplitude:" << ampl << " ...\n";
    THROW(fatal_error, "Process not found");
  }
  pit->second->Differential(*campl, 2 | 4 | 128);
  campl->Delete();
  std::vector<std::vector<Complex>> cols;
  pit->second->FillAmplitudes(spinampls, cols);
}

void Comix_Interface::InitializeProcesses(const Cluster_Amplitude_PM& ampls)
{
  DEBUG_FUNC("");
  auto& s = Settings::GetMainSettings();
  const auto graph_path =
      s["PRINT_EWSUDAKOV_GRAPHS"].SetDefault("").Get<std::string>();
  for (const auto& kv : ampls) {
    const auto& ampl = kv.second;
    msg_Debugging() << "Initialize process for ampl=" << *ampl << std::endl;
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
  PHASIC::Process_Base* proc =
      p_proc->Generator()->Generators()->InitializeProcess(pi, false);
  if (proc == NULL){
    msg_Error() << "Invalid process: " << pi << std::endl;
    THROW(fatal_error, "Invalid process");
  }
  proc->SetSelector(Selector_Key{});
  proc->SetScale(Scale_Setter_Arguments(
      MODEL::s_model, "VAR{" + ToString(sqr(rpa->gen.Ecms())) + "}",
      "Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  proc->Get<COMIX::Process_Base>()->Tests();
  proc->FillProcessMap(&m_apmap);
}
