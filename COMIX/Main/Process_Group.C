#include "COMIX/Main/Process_Group.H"

#include "COMIX/Main/Single_Process.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

COMIX::Process_Group::Process_Group(): 
  COMIX::Process_Base(this) {}

COMIX::Process_Group::Process_Group(Model *const model):
  COMIX::Process_Base(this,model) {}

bool COMIX::Process_Group::Initialize(std::map<std::string,std::string> *const pmap,
				      std::vector<Single_Process*> *const procs)
{
  return Process_Base::Initialize(pmap,procs);
}

PHASIC::Process_Base *COMIX::Process_Group::GetProcess(const PHASIC::Process_Info &pi) const
{
  return new Single_Process();
}

bool COMIX::Process_Group::Initialize(PHASIC::Process_Base *const proc)
{
  COMIX::Process_Base *cdxs(proc->Get<COMIX::Process_Base>());
  cdxs->SetModel(p_model);
  cdxs->SetPSMC(m_psmc);
  cdxs->SetGPath(m_gpath);
  proc->Integrator()->SetHelicityScheme(p_int->HelicityScheme());
  if (!cdxs->Initialize(p_pmap,p_umprocs)) return false;
  if (!cdxs->MapProcess())
    if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  proc->SetParent((PHASIC::Process_Base*)this);
  return true;
}

bool COMIX::Process_Group::MapProcess()
{
  return false;
}

bool COMIX::Process_Group::GeneratePoint()
{
  bool zero=true;
  m_last=0.0;
  for (size_t i(0);i<m_procs.size();++i)
    if (m_procs[i]->GeneratePoint()) zero=false;
  return !zero;
}

bool COMIX::Process_Group::Tests()
{
  for (size_t i=0;i<m_procs.size();++i)
    if (!m_procs[i]->Get<COMIX::Process_Base>()->Tests()) return false;
  return true;
}

void COMIX::Process_Group::InitPSGenerator(const size_t &ismode)
{
  if (!(ismode&1)) {
    p_psgen = new PS_Generator(this);
  }
  else {
    for (size_t i(0);i<Size();++i) 
      (*this)[i]->Get<COMIX::Process_Base>()->InitPSGenerator(ismode);
  }
}

void COMIX::Process_Group::ConstructPSVertices(PS_Generator *ps)
{
  for (size_t i(0);i<m_procs.size();++i)
    m_procs[i]->Get<COMIX::Process_Base>()->ConstructPSVertices(ps);
}

bool COMIX::Process_Group::FillIntegrator(Phase_Space_Handler *const psh)
{
  return Process_Base::FillIntegrator(psh);
}

void COMIX::Process_Group::UpdateIntegrator(Phase_Space_Handler *const psh)
{
  Process_Base::UpdateIntegrator(psh);
}
