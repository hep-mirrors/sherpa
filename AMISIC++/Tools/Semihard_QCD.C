#include "AMISIC++/Tools/Semihard_QCD.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/ISR_Vegas.H"
#include "PDF/Main/ISR_Handler.H"

using namespace AMISIC;

Semihard_QCD::Semihard_QCD(ATOOLS::Default_Reader *const read):
  ME_Generator_Base("Amisic")
{
  SetPSMasses(read);
  SetFSRInterface(NULL);
  SetFSRMode(0);
}

Semihard_QCD::~Semihard_QCD()
{
}

void Semihard_QCD::CreateFSRChannels() 
{
  if (m_fsrmode==0 || p_fsrinterface==NULL) {
    p_int->PSHandler()->FSRIntegrator()->DropAllChannels();
    p_int->PSHandler()->FSRIntegrator()->
      Add(new PHASIC::S1Channel(m_nin,m_nout,&m_flavs.front()));
    p_int->PSHandler()->FSRIntegrator()->
      Add(new PHASIC::T1Channel(m_nin,m_nout,&m_flavs.front()));
    p_int->PSHandler()->FSRIntegrator()->
      Add(new PHASIC::U1Channel(m_nin,m_nout,&m_flavs.front()));
    m_fsrmode=1;
  }
  else {
    if (m_fsrmode==3) {
      p_int->PSHandler()->FSRIntegrator()->DropAllChannels(false);
      m_fsrmode=0;
    }
    if (m_fsrmode==2) {
      p_int->PSHandler()->FSRIntegrator()->DropAllChannels();
      p_int->PSHandler()->FSRIntegrator()->Add(p_fsrinterface);
      p_fsrinterface->SetAlpha(1.0);
      m_fsrmode=1;
    }
  }
}

void Semihard_QCD::CreateISRChannels() 
{
  PHASIC::Multi_Channel *isr=p_int->PSHandler()->ISRIntegrator();
  isr->DropAllChannels();
  PHASIC::Single_Channel *channel =
    new PHASIC::Simple_Pole_Uniform_V(1.0," isr",p_int->PSHandler()->GetInfo());
  channel->SetAlpha(1.0);
  isr->Add(channel);
  isr->Reset();
}

void Semihard_QCD::InitIntegrators() 
{
  p_int->PSHandler()->CreateIntegrators();
  p_int->PSHandler()->ISRIntegrator()->Reset();
  p_int->PSHandler()->FSRIntegrator()->Reset();
}
      
