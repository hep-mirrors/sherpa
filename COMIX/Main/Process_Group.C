#include "COMIX/Main/Process_Group.H"

#include "COMIX/Main/Single_Process.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "COMIX/Phasespace/PS_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/MathTools.H"
#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"

using namespace COMIX;
using namespace PHASIC;
using namespace ATOOLS;

COMIX::Process_Group::Process_Group(): 
  COMIX::Process_Base(this), m_nproc(0) {}

COMIX::Process_Group::Process_Group(MODEL::Model_Base *const model):
  COMIX::Process_Base(this,model), m_nproc(0) {}

bool COMIX::Process_Group::Initialize(std::map<std::string,std::string> *const pmap,
				      std::vector<Single_Process*> *const procs)
{
  std::string mapfile(rpa->gen.Variable("SHERPA_CPP_PATH")
		      +"/Process/Comix/"+m_name+".map");
  msg_Debugging()<<"checking for '"<<mapfile<<"' ... "<<std::flush;
  if (!FileExists(mapfile,1)) {
    msg_Debugging()<<"not found"<<std::endl;
  }
  else {
    msg_Debugging()<<"found"<<std::endl;
    My_In_File map(mapfile);
    if (map.Open()) {
      while (!map->eof()) {
	std::string src, dest;
	*map>>src>>dest;
	if (src!="" && dest!="x")
	  THROW(fatal_error,"Corrupted map file '"+mapfile+"'");
	if (map->eof()) break;
	(*pmap)[src]=dest;
	msg_Debugging()<<" map '"<<src<<"' onto '"<<dest<<"'\n";
      }
    }
  }
  return COMIX::Process_Base::Initialize(pmap,procs,m_blocks,m_nproc);
}

PHASIC::Process_Base *COMIX::Process_Group::GetProcess(const PHASIC::Process_Info &pi) const
{
  return new Single_Process();
}

bool COMIX::Process_Group::Initialize(PHASIC::Process_Base *const proc)
{
  COMIX::Process_Base *cdxs(proc->Get<COMIX::Process_Base>());
  cdxs->SetModel(p_model);
  cdxs->SetGPath(m_gpath);
  cdxs->SetCTS(p_cts);
  proc->Integrator()->SetHelicityScheme(p_int->HelicityScheme());
  proc->SetParent((PHASIC::Process_Base*)this);
  if (!cdxs->Initialize(p_pmap,p_umprocs,m_blocks,m_nproc)) return false;
  if (s_partcommit) My_In_File::CloseDB
    (rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/",0);
  if (!cdxs->MapProcess())
    if (!msg_LevelIsTracking()) msg_Info()<<"."<<std::flush;
  return true;
}

bool COMIX::Process_Group::MapProcess()
{
  return false;
}

bool COMIX::Process_Group::GeneratePoint()
{
  bool zero=true;
  m_last["Main"]=0.0;
  m_last["All"]=0.0;
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
    p_psgen = std::make_shared<PS_Generator>(this);
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
  bool res(COMIX::Process_Base::FillIntegrator(psh));
  for (size_t i(0);i<m_procs.size();++i) {
    m_procs[i]->Get<COMIX::Process_Base>()->SetISMC(p_ismc);
    m_procs[i]->Get<COMIX::Process_Base>()->SetFSMC(p_fsmc);
  }
  if ((m_pinfo.m_fi.NLOType()&nlo_type::vsub) && p_ismc) {
    p_ismc->AddERan("z_1");
    p_ismc->AddERan("z_2");
  }
  return res;
}
