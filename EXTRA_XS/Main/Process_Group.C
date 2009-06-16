#include "EXTRA_XS/Main/Process_Group.H"
#include "EXTRA_XS/Main/Single_Process.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/FSR_Channel.H"

using namespace ATOOLS;
using namespace EXTRAXS;
using namespace PHASIC;

EXTRAXS::Process_Group::Process_Group()
{
}

EXTRAXS::Process_Group::~Process_Group()
{
}

PHASIC::Process_Base *EXTRAXS::Process_Group::GetProcess(const PHASIC::Process_Info &pi) const
{
  return new Single_Process();
}

bool EXTRAXS::Process_Group::Initialize(Process_Base *const proc)
{
  if (!proc->Get<Single_Process>()->Initialize()) return false;
  proc->SetParent((Process_Base*)this);
  if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  return true;
}

bool EXTRAXS::Process_Group::FillIntegrator
(PHASIC::Phase_Space_Handler *const psh)
{
  Multi_Channel *mc(psh->FSRIntegrator());
  mc->DropAllChannels();
  size_t sintt(0);
  for (size_t i(0);i<m_procs.size();++i)
    if (m_procs[i]->Get<EXTRAXS::Single_Process>()->GetME()) 
      sintt|=m_procs[i]->Get<EXTRAXS::Single_Process>()->GetME()->SIntType();
    else sintt|=7;
  if (sintt&1) mc->Add(new S1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  if (sintt&2) mc->Add(new T1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  if (sintt&4) mc->Add(new U1Channel(m_nin,m_nout,(Flavour*)&Flavours().front()));
  return false;
}      
