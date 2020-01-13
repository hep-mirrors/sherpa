#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "ATOOLS/Phys/Color.H"
#include "COMIX/Main/Single_Process.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"
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
                                 EWSudakov_Amplitudes& ampls)
    : p_proc{proc}, p_model_he{nullptr}, p_model{nullptr}
{
  InitializeHighEnergyModel();
  InitializeProcesses(ampls);
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
  if (pit->second == NULL)
    THROW(fatal_error, "Process not found");
  double me2 = pit->second->Differential(*campl, 2 | 4 | 128);
  PRINT_VAR(me2);
  campl->Delete();
  std::vector<std::vector<Complex>> cols;
  pit->second->FillAmplitudes(spinampls, cols);

  const auto loprocmapit_he = m_apmap_he.find(nlo_type::lo);
  if (loprocmapit_he == m_apmap_he.end())
    THROW(fatal_error, "LO entry in process map not found");
  Cluster_Amplitude* campl_he(ampl.Copy());
  campl_he->SetMuR2(sqr(rpa->gen.Ecms()));
  campl_he->SetMuF2(sqr(rpa->gen.Ecms()));
  campl_he->SetMuQ2(sqr(rpa->gen.Ecms()));
  pname = Process_Base::GenerateName(campl_he);
  PRINT_VAR(pname);
  pit = loprocmapit_he->second->find(pname);
  if (pit->second == NULL)
    THROW(fatal_error, "Process not found");
  me2 = pit->second->Differential(*campl_he, 2 | 4 | 128);
  PRINT_VAR(me2);
  campl_he->Delete();
}

void Comix_Interface::InitializeProcesses(EWSudakov_Amplitudes& ampls)
{
  p_model = s_model;
  DEBUG_FUNC("");
  auto& s = Settings::GetMainSettings();
  const auto graph_path =
      s["PRINT_EWSUDAKOV_GRAPHS"].SetDefault("").Get<std::string>();
  for (auto& kv : ampls) {
    const auto& ampl = kv.second;
    msg_Debugging() << "Initialize process for ampl=" << *ampl << std::endl;
    const Process_Info pi = CreateProcessInfo(ampl, graph_path);
    InitializeProcess(pi);
  }

  MODEL::Model_Base* model = s_model;
  MODEL::s_model = p_model_he;

  p_model = p_model_he;
  for (auto& kv : ampls) {
    const auto& ampl = kv.second;
    msg_Debugging() << "Initialize HE process for ampl=" << *ampl << std::endl;
    const Process_Info pi = CreateProcessInfo(ampl, graph_path, "Sudakov_HE");
    InitializeProcess(pi);
  }
  p_model = model;

  MODEL::s_model = model;
}

Process_Info
Comix_Interface::CreateProcessInfo(const Cluster_Amplitude_UP& ampl,
                                   const std ::string& graph_path,
                                   const std::string& suffix)
{
  Process_Info pi;
  pi.m_addname = "__" + suffix;
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
  p_proc->Generator()->Generators()->SetModel(p_model);
  PHASIC::Process_Base* proc =
      p_proc->Generator()->Generators()->InitializeProcess(pi, false);
  if (proc == NULL)
    THROW(fatal_error, "Invalid process");
  proc->SetSelector(Selector_Key{});
  proc->SetScale(Scale_Setter_Arguments(
      MODEL::s_model, "VAR{" + ToString(sqr(rpa->gen.Ecms())) + "}",
      "Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  proc->Get<COMIX::Process_Base>()->Tests();
  proc->FillProcessMap((pi.m_addname == "__Sudakov_HE") ? &m_apmap_he : &m_apmap);
}

bool Comix_Interface::InitializeHighEnergyModel()
{
  Settings& s = Settings::GetMainSettings();
  s["SIN2THETAW"].OverrideScalar(0.3);
  if (p_model_he) delete p_model_he;
  std::string name(s["MODEL"].Get<std::string>());
  p_model_he=Model_Base::Model_Getter_Function::
    GetObject(name, Model_Arguments(true));
  //if (p_model_he==NULL) {
  //  if (!s_loader->LoadLibrary("Sherpa"+name))
  //    THROW(missing_module,"Cannot load model library Sherpa"+name+".");
  //  p_model_he=Model_Base::Model_Getter_Function::
  //    GetObject(name, Model_Arguments(true));
  //}
  if (p_model_he==NULL) THROW(not_implemented,"Model not implemented");

  PDF::ISR_Handler_Map isrmap;
  isrmap[PDF::isr::id::hard_subprocess] = p_proc->Integrator()->ISR();
  if (!p_model_he->ModelInit(*SHERPA::s_inithandler->GetISRHandlers()))
    THROW(critical_error,"Model cannot be initialized");
  p_model_he->InitializeInteractionModel();
  return 1;
}

