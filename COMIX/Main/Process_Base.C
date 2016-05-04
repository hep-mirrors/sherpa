#include "COMIX/Main/Process_Base.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "PHASIC++/Channels/VHAAG.H"
#include "COMIX/Phasespace/PS_Channel.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

std::string COMIX::ComixLogo()
{
  if (!msg->Modifiable()) return "Comix";
  return "\033[31mC\033[32mo\033[34mm\033[0mi\033[33mx\033[0m";
}

int COMIX::Process_Base::s_partcommit=0;

COMIX::Process_Base::Process_Base
(PHASIC::Process_Base *const proc,MODEL::Model_Base *const model):
  p_proc(proc), p_model(model), p_psgen(NULL),
  m_cls(-1), m_hls(-1), p_cts(NULL),
  p_pmap(NULL), p_umprocs(NULL) {}

COMIX::Process_Base::~Process_Base() 
{
}

bool COMIX::Process_Base::Initialize(std::map<std::string,std::string> *const pmap,
				     std::vector<Single_Process*> *const procs)
{
  p_pmap=pmap;
  p_umprocs=procs;
  p_proc->Integrator()->SetColorScheme(cls::sample);
  return true;
}

bool COMIX::Process_Base::FillIntegrator(Phase_Space_Handler *const psh)
{
  if (p_proc->NOut()==1) return false;
  Multi_Channel *mc(psh->FSRIntegrator());
  mc->DropAllChannels();
  PS_Channel *ch(new PS_Channel(p_proc->NIn(),p_proc->NOut(),
				(Flavour*)&p_proc->Flavours().front(),this));
  InitPSGenerator(0);
  mc->Add(ch);
  return false;
}      
