#include "PHASIC++/EWSudakov/HE_Comix_Interface.H"

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

std::unique_ptr<MODEL::Model_Base> HE_Comix_Interface::p_model_he {nullptr};

HE_Comix_Interface::HE_Comix_Interface(Process_Base* proc,
                                       EWSudakov_Amplitudes& ampls)
    : Comix_Interface {proc, "Sudakov_HE"}
{
  InitializeHighEnergyModel();
  InitializeProcesses(ampls);
}

void HE_Comix_Interface::InitializeProcesses(EWSudakov_Amplitudes& ampls)
{
  Model_Base* model = s_model;
  s_model = p_model_he.get();
  Comix_Interface::InitializeProcesses(ampls);
  s_model = model;
}

void HE_Comix_Interface::InitializeProcess(const Process_Info& pi)
{
  p_proc->Generator()->Generators()->SetModel(p_model_he.get());
  Comix_Interface::InitializeProcess(pi);
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

  // TODO: probably we want to suppress this output at some point, since most
  // (all?) of it is just duplicating the "normal" model init output
  msg_Out() << '\n';
  PRINT_FUNC("");

  // create model
  Settings& s = Settings::GetMainSettings();
  std::string name(s["MODEL"].Get<std::string>());
  p_model_he.reset(Model_Base::Model_Getter_Function::GetObject(
      name, Model_Arguments(true)));

  // TODO: re-consider if this needs to be re-enabled; but probably not, since
  // we only ever use SMGold in conjunction with the EW Sudakovs; but then this
  // should be enforced and made more explicit here
  //if (p_model_he==NULL) {
  //  if (!s_loader->LoadLibrary("Sherpa"+name))
  //    THROW(missing_module,"Cannot load model library Sherpa"+name+".");
  //  p_model_he=Model_Base::Model_Getter_Function::
  //    GetObject(name, Model_Arguments(true));
  //}

  if (p_model_he == nullptr)
    THROW(not_implemented, "Model not implemented");

  // init model
  // TODO: devise another way to get the ISR handlers and undo making the
  // Initialization_Handler a global object; the best way might be to store the
  // necessary information to create a copy within model, including the ISR
  // handlers used to initialize it
  if (!p_model_he->ModelInit(*SHERPA::s_inithandler->GetISRHandlers()))
    THROW(critical_error, "Model cannot be initialized");
  p_model_he->InitializeInteractionModel();
  return true;
}
