#include "AddOns/EWSud/HE_Comix_Interface.H"

#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "SHERPA/Main/Sherpa.H"

using namespace ATOOLS;
using namespace MODEL;
using namespace PHASIC;
using namespace EWSud;

std::unique_ptr<MODEL::Model_Base> HE_Comix_Interface::p_model_he {nullptr};
NLOTypeStringProcessMap_Map HE_Comix_Interface::s_apmap_he;

HE_Comix_Interface::HE_Comix_Interface(Process_Base* proc,
                                       const Amplitudes& ampls)
    : Comix_Interface {proc, "Sudakov_HE"}
{
  InitializeHighEnergyModel();
  InitializeProcesses(ampls.GoldstoneOnly());
}

void HE_Comix_Interface::InitializeProcesses(const Cluster_Amplitude_PM& ampls)
{
  Model_Base* model = s_model;
  s_model = p_model_he.get();
  p_proc->Generator()->Generators()->SetModel(p_model_he.get());
  Comix_Interface::InitializeProcesses(ampls);
  s_model = model;
  p_proc->Generator()->Generators()->SetModel(s_model);
}

void HE_Comix_Interface::ResetWithEWParameters(const MODEL::EWParameters& p)
{
  p_model_he->ResetVerticesWithEWParameters(p);
}

bool HE_Comix_Interface::InitializeHighEnergyModel()
{
  static bool did_initialize {false};
  if (did_initialize) {
    return true;
  } else {
    did_initialize = true;
  }

  // Suppress model initialization output, since most
  // (all?) of it is just duplicating the "normal" model init output
  auto level = msg->Level();
  msg->SetLevel(0);

  // create model
  Settings& s = Settings::GetMainSettings();
  std::string name(s["MODEL"].Get<std::string>());
  p_model_he.reset(Model_Base::Model_Getter_Function::GetObject(
      name, Model_Arguments(true)));
  if (!p_model_he)
    THROW(missing_module,"Cannot load model library Sherpa"+name+".");

  // init model
  if (!p_model_he->ModelInit(s_model->ISRHandlerMap()))
    THROW(critical_error, "Model cannot be initialized");
  p_model_he->InitializeInteractionModel();

  msg->SetLevel(level);

  return true;
}

NLOTypeStringProcessMap_Map& HE_Comix_Interface::ProcessMap()
{
  return s_apmap_he;
}

const NLOTypeStringProcessMap_Map& HE_Comix_Interface::ProcessMap() const
{
  return s_apmap_he;
}
