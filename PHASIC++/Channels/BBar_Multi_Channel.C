#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/ISR_Handler.H"

using namespace PHASIC;
using namespace ATOOLS;

BBar_Multi_Channel::BBar_Multi_Channel
(Process_Base *const proc,Process_Base *const sproc,
 Phase_Space_Handler *const psh):
  Multi_Channel("BBar_MC"), p_proc(proc),
  p_fsmc(psh->FSRIntegrator()),
  p_cuts(p_proc->Integrator()->PSHandler()->Cuts())
{
  DEBUG_FUNC(p_proc->Name());
  nin=p_proc->NIn();
  nout=p_proc->NOut();
  m_eeg.InitDipoles(p_proc,sproc,psh);
}

BBar_Multi_Channel::~BBar_Multi_Channel()
{
  delete p_fsmc;
}

Dipole_Params BBar_Multi_Channel::Active(Process_Base *const bviproc) const
{
  return m_eeg.Active(bviproc);
}

void BBar_Multi_Channel::Reset() 
{
  p_fsmc->Reset();
  Print();
}

void BBar_Multi_Channel::GenerateWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  m_eeg.GenerateWeight(p_cuts,true);
  p_fsmc->GenerateWeight(p,cuts);
  m_weight=p_fsmc->Weight();
}

void BBar_Multi_Channel::GeneratePoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  p_fsmc->GeneratePoint(p,cuts);
  m_eeg.GeneratePoint(Vec4D_Vector(p,&p[nin+nout]),cuts);
}

void BBar_Multi_Channel::GenerateEmissionPoint
(const ATOOLS::Cluster_Amplitude &ampl)
{
  Vec4D_Vector p(nin+nout);
  for (size_t i(0);i<nin+nout;++i)
    p[i]=i<nin?-ampl.Leg(i)->Mom():ampl.Leg(i)->Mom();
  m_eeg.GeneratePoint(p,p_cuts);
  m_eeg.GenerateWeight(p_cuts,true);
}

void BBar_Multi_Channel::AddPoint(double value)
{ 
  m_lastdice=-1;
  Multi_Channel::AddPoint(value);
  p_fsmc->AddPoint(value*p_proc->LastB()/p_proc->Last());
  m_eeg.AddPoint(value/p_proc->Last());
}

void BBar_Multi_Channel::Optimize(double error)
{ 
  p_fsmc->Optimize(error);
  m_eeg.Optimize();
  Print();
}

void BBar_Multi_Channel::EndOptimize(double error)
{ 
  p_fsmc->EndOptimize(error);
  m_eeg.EndOptimize();
}

void BBar_Multi_Channel::MPISync()
{
  Multi_Channel::MPISync();
  p_fsmc->MPISync();
  m_eeg.MPISync();
}

bool BBar_Multi_Channel::OptimizationFinished()
{
  return p_fsmc->OptimizationFinished();
}

void BBar_Multi_Channel::WriteOut(std::string pid)
{ 
  Multi_Channel::WriteOut(pid+"_BBMC");
  p_fsmc->WriteOut(pid);
  m_eeg.WriteOut(pid);
}
    
bool BBar_Multi_Channel::ReadIn(std::string pid)
{
  Multi_Channel::ReadIn(pid+"_BBMC");
  if (!p_fsmc->ReadIn(pid)) return false;
  if (!m_eeg.ReadIn(pid)) return false;
  return true;
}

size_t BBar_Multi_Channel::Number()
{
  return p_fsmc->Number();
}

void BBar_Multi_Channel::ISRInfo
(int i,int &t,double &m,double &w)
{
  p_fsmc->ISRInfo(i,t,m,w);
}

void BBar_Multi_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,
 std::vector<double> &ws) const
{
  p_fsmc->ISRInfo(ts,ms,ws);
}

std::string BBar_Multi_Channel::Name()
{
  return name;
}

std::string BBar_Multi_Channel::ChID()
{
  return name;
}

void BBar_Multi_Channel::Print() 
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
