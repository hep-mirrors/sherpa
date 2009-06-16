#include "PHASIC++/Process/Single_Process.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

Single_Process::Single_Process()
{
}

Single_Process::~Single_Process()
{
}

size_t Single_Process::Size() const
{
  return 1;
}

Process_Base *Single_Process::operator[](const size_t &i)
{
  if (i==0) return this;
  return NULL;
}

void Single_Process::DeSelect()
{
  p_selected=NULL;
}

bool Single_Process::SelectOne()
{
  p_selected=this;
  return true;
}

Weight_Info *Single_Process::OneEvent() 
{
  SelectOne();
  return p_int->PSHandler()->OneEvent(this);
}

Weight_Info *Single_Process::WeightedEvent(const int mode) 
{
  SelectOne();
  return p_int->PSHandler()->WeightedEvent(this,mode);
}

double Single_Process::KFactor()
{
  return p_kfactor->KFactor();
}

double Single_Process::KFactor2()
{
  return p_kfactor->KFactor2(p_scale->Scale2());
}

bool Single_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
          m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
        p_int->ISR()->SetPartonMasses(&m_flavs.front());
      }
    }
  }
  psh->InitCuts();
  if (p_int->ISR())
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"
            <<m_name<<"' ("<<(p_gen?p_gen->Name():"")<<")"<<std::endl;
  double totalxs(psh->Integrate()/rpa.Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->TotalXS()>0.0) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
  }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}

void Single_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
}
