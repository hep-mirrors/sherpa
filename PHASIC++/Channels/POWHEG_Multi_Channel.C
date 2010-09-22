#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/POWHEG_Multi_Channel.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/ISR_Handler.H"

using namespace PHASIC;
using namespace ATOOLS;

POWHEG_Multi_Channel::POWHEG_Multi_Channel
(POWHEG_Process *const proc,Process_Base *const sproc,
 Phase_Space_Handler *const psh):
  Multi_Channel("POWHEG_MC"), p_proc(proc),
  p_fsmc(psh->FSRIntegrator()), p_ismc(psh->ISRIntegrator()),
  p_bampl(Cluster_Amplitude::New()),
  p_cuts(proc->Integrator()->PSHandler()->Cuts()),
  m_ismode(psh->Process()->ISR()->On()),
  m_emode(proc->Info().Has(nlo_type::real))
{
  nin=p_proc->NIn();
  nout=p_proc->NOut();
  p_bampl->SetNIn(nin);
  for (int i(0);i<nin+nout;++i)
    p_bampl->CreateLeg(Vec4D(),Flavour(kf_jet));
  Vec4D_Vector p(nin+nout);
  do {
    psh->TestPoint(&p.front(),&sproc->Info());
    for (int i(0);i<nin+nout;++i)
      p_bampl->Leg(i)->SetMom(i<nin?-p[i]:p[i]);
    sproc->Differential(*p_bampl);
  } while (!p_proc->InitSubtermInfo());
  p_bampl->Legs().back()->Delete();
  p_bampl->Legs().pop_back();
  m_eeg.InitDipoles(p_proc,sproc,psh);
  m_pb.resize(nin+nout-1);
  m_pr.resize(nin+nout);
  m_isrspkey.Assign("s' isr",4,0,psh->GetInfo());
  m_isrykey.Assign("y isr",3,0,psh->GetInfo());
}

POWHEG_Multi_Channel::~POWHEG_Multi_Channel()
{
  p_bampl->Delete();
  delete p_fsmc;
}

void POWHEG_Multi_Channel::Reset() 
{
  p_fsmc->Reset();
  Print();
}

double POWHEG_Multi_Channel::GenerateWeight
(const Cluster_Amplitude &ampl) 
{
  Vec4D_Vector p(ampl.Legs().size());
  for (size_t i(0);i<ampl.NIn();++i) p[i]=-ampl.Leg(i)->Mom();
  for (size_t i(ampl.NIn());i<p.size();++i) p[i]=ampl.Leg(i)->Mom();
  if (p.size()==size_t(nin+nout)) {
    m_eeg.GenerateWeight(p,p_cuts);
    return m_eeg.Weight();
  }
  p_fsmc->GenerateWeight(&p.front(),p_cuts);
  if (p_ismc==NULL) return p_fsmc->Weight();
  m_isrspkey[3]=(m_pb[0]+m_pb[1]).Abs2();
  m_isrykey[2]=(m_pb[0]+m_pb[1]).Y();
  p_ismc->GenerateWeight(m_ismode);
  return p_ismc->Weight()*p_fsmc->Weight();
}

void POWHEG_Multi_Channel::GenerateWeight
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  if (cuts==NULL) {
    m_eeg.GenerateWeight(m_pr,p_cuts,true);
  }
  else {
    if (m_emode) m_eeg.GenerateWeight(m_pr,cuts);
    else p_fsmc->GenerateWeight(&m_pb.front(),cuts);
  }
  m_bweight=p_fsmc->Weight();
  m_weight=m_emode?m_eeg.Weight():m_bweight;
  if (p_ismc && m_emode) {
    m_weight/=p_ismc->Weight();
  }
}

void POWHEG_Multi_Channel::GeneratePoint
(ATOOLS::Vec4D *p,Cut_Data *cuts)
{
  p_fsmc->GeneratePoint(p,cuts);
  for (int i(0);i<nin+nout-1;++i) {
    m_pb[i]=p[i];
    p_bampl->Leg(i)->SetMom(i<nin?-p[i]:p[i]);
  }
  m_pr=m_eeg.GeneratePoint(m_pb,cuts);
  for (int i(0);i<nin+nout;++i) p[i]=m_pr[i];
}

void POWHEG_Multi_Channel::GenerateEmissionPoint
(const Vec4D_Vector& pb)
{
  for (int i(0);i<nin+nout-1;++i) {
    m_pb[i]=pb[i];
    p_bampl->Leg(i)->SetMom(i<nin?-pb[i]:pb[i]);
  }
  m_pr=m_eeg.GeneratePoint(m_pb,p_cuts);
}

void POWHEG_Multi_Channel::AddPoint(double value)
{ 
  m_lastdice=-1;
  Multi_Channel::AddPoint(value);
  p_fsmc->AddPoint(value/p_proc->Last()*
		   (p_proc->LastB()+p_proc->LastVI()+p_proc->LastRS()));
  if (m_emode) m_eeg.AddPoint(value/p_proc->Last()*p_proc->LastRS());
}

void POWHEG_Multi_Channel::Optimize(double error)
{ 
  p_fsmc->Optimize(error);
  if (m_emode) m_eeg.Optimize();
  Print();
}

void POWHEG_Multi_Channel::EndOptimize(double error)
{ 
  p_fsmc->EndOptimize(error);
  if (m_emode) m_eeg.EndOptimize();
}

bool POWHEG_Multi_Channel::OptimizationFinished()
{
  return p_fsmc->OptimizationFinished();
}

void POWHEG_Multi_Channel::WriteOut(std::string pid)
{ 
  p_fsmc->WriteOut(pid);
  if (m_emode) m_eeg.WriteOut(pid);
}
    
bool POWHEG_Multi_Channel::ReadIn(std::string pid)
{
  if (!p_fsmc->ReadIn(pid)) return false;
  if (m_emode) m_eeg.ReadIn(pid);
  return true;
}

size_t POWHEG_Multi_Channel::Number()
{
  return p_fsmc->Number();
}

void POWHEG_Multi_Channel::ISRInfo
(int i,int &t,double &m,double &w)
{
  p_fsmc->ISRInfo(i,t,m,w);
}

void POWHEG_Multi_Channel::ISRInfo
(std::vector<int> &ts,std::vector<double> &ms,
 std::vector<double> &ws) const
{
  p_fsmc->ISRInfo(ts,ms,ws);
}

std::string POWHEG_Multi_Channel::Name()
{
  return name;
}

std::string POWHEG_Multi_Channel::ChID()
{
  return name;
}

size_t POWHEG_Multi_Channel::ActiveIdx() const
{
  return m_eeg.Active()->Idx();
}

double POWHEG_Multi_Channel::SelectionWeight(const size_t &idx) const
{
  return m_eeg.SelectionWeight(idx);
}

void POWHEG_Multi_Channel::Print() 
{
  if (!msg_LevelIsTracking()) return;
  msg_Out()<<name<<" {\n";
  {
    msg_Indent();
    p_fsmc->Print();
    if (m_emode) m_eeg.Print();
  }
  msg_Out()<<"}\n";
}                 
