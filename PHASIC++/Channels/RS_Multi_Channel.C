#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/RS_Multi_Channel.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace PHASIC;
using namespace ATOOLS;

RS_Multi_Channel::RS_Multi_Channel
(Process_Base *const proc,Phase_Space_Handler *const psh):
  Multi_Channel("RS_MC"), p_proc(proc),
  p_fsmc(psh->FSRIntegrator())
{
  DEBUG_FUNC(p_proc->Name());
  nin=p_proc->NIn();
  nout=p_proc->NOut();
  m_eeg.InitDipoles(proc,psh);
}

RS_Multi_Channel::~RS_Multi_Channel()
{
  delete p_fsmc;
}

void RS_Multi_Channel::Reset() 
{
  Multi_Channel::Reset();
  p_fsmc->Reset();
  Print();
}

bool RS_Multi_Channel::Initialize() 
{
  DEBUG_FUNC("");
  p_fsmc->SetNout(nout-1);
  bool res(p_fsmc->Initialize());
  m_incisr=p_fsmc->IncludesISR();
  return res;
}

void RS_Multi_Channel::Add(Single_Channel *channel) 
{
  return p_fsmc->Add(channel);
}

void RS_Multi_Channel::GenerateWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts,bool compute)
{
  if (!m_eeg.GenerateWeight(cuts,true)) m_weight=0;
  else m_weight=m_eeg.Weight();
}

void RS_Multi_Channel::GeneratePoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  p_fsmc->GeneratePoint(p,cuts);
  m_status=p_fsmc->Status();
  if (!m_status) return;
  m_status&=m_eeg.GeneratePoint(Vec4D_Vector(p,&p[nin+nout-1]),cuts);
  if (!m_status) return;
  for (size_t i(0);i<nin+nout;++i) p[i]=m_eeg.Momenta()[i];
  m_sprime=(p[0]+p[1]).Abs2();
  m_y=(p[0]+p[1]).Y();
}

void RS_Multi_Channel::AddPoint(double value)
{ 
  m_lastdice=-1;
  Multi_Channel::AddPoint(value);
  p_fsmc->AddPoint(value);
  m_eeg.AddPoint(value);
}

void RS_Multi_Channel::Optimize(double error)
{ 
  p_fsmc->Optimize(error);
  m_eeg.Optimize();
  Print();
}

void RS_Multi_Channel::EndOptimize(double error)
{ 
  p_fsmc->EndOptimize(error);
  m_eeg.EndOptimize();
}

void RS_Multi_Channel::MPISync()
{
  Multi_Channel::MPISync();
  p_fsmc->MPISync();
  m_eeg.MPISync();
}

bool RS_Multi_Channel::OptimizationFinished()
{
  return p_fsmc->OptimizationFinished();
}

void RS_Multi_Channel::WriteOut(std::string pid)
{ 
  Multi_Channel::WriteOut(pid+"_BBMC");
  p_fsmc->WriteOut(pid);
  m_eeg.WriteOut(pid);
}
    
bool RS_Multi_Channel::ReadIn(std::string pid)
{
  Multi_Channel::ReadIn(pid+"_BBMC");
  if (!p_fsmc->ReadIn(pid)) return false;
  if (!m_eeg.ReadIn(pid)) return false;
  return true;
}

size_t RS_Multi_Channel::Number()
{
  return p_fsmc->Number();
}

void RS_Multi_Channel::ISRInfo
(int i,int &t,double &m,double &w)
{
  p_fsmc->ISRInfo(i,t,m,w);
}

void RS_Multi_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,
 std::vector<double> &ws) const
{
  p_fsmc->ISRInfo(ts,ms,ws);
}

std::string RS_Multi_Channel::Name()
{
  return name;
}

std::string RS_Multi_Channel::ChID()
{
  return name;
}

void RS_Multi_Channel::Print() 
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<name<<" {\n";
  {
    msg_Indent();
    p_fsmc->Print();
    m_eeg.Print();
  }
  msg_Out()<<"}\n";
}                 
