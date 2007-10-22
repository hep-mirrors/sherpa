#include "Semihard_QCD.H"

#include "Phase_Space_Handler.H"
#include "Regulator_Base.H"
#include "FSR_Channel.H"
#include "ISR_Vegas.H"
#include "ISR_Handler.H"

#define INIFL ATOOLS::Flavour(ATOOLS::kf::jet)

using namespace AMISIC;

Semihard_QCD::Semihard_QCD(BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr,
			   ATOOLS::Selector_Data *const seldata,
			   const ATOOLS::Flavour *flavours,
			   const int scalescheme,
			   const int kfactorscheme):
  XS_Group(2,2,flavours,(PHASIC::scl::scheme)(scalescheme),kfactorscheme,beam,isr,seldata,NULL)
{
  SetFSRInterface(NULL);
  SetFSRMode(0);
}

Semihard_QCD::~Semihard_QCD()
{
}

void Semihard_QCD::CreateFSRChannels() 
{
  if (m_fsrmode==0 || p_fsrinterface==NULL) {
    p_pshandler->FSRIntegrator()->DropAllChannels();
    if (p_isrhandler->KMROn()>0) {
      p_pshandler->FSRIntegrator()->
	Add(new PHASIC::T2Channel(m_nin,m_nout,p_flavours));
      p_pshandler->FSRIntegrator()->
	Add(new PHASIC::T3Channel(m_nin,m_nout,p_flavours));
    }
    else {
      p_pshandler->FSRIntegrator()->
	Add(new PHASIC::S1Channel(m_nin,m_nout,p_flavours));
      p_pshandler->FSRIntegrator()->
	Add(new PHASIC::T1Channel(m_nin,m_nout,p_flavours));
      p_pshandler->FSRIntegrator()->
	Add(new PHASIC::U1Channel(m_nin,m_nout,p_flavours));
    }
    m_fsrmode=1;
  }
  else {
    if (m_fsrmode==3) {
      p_pshandler->FSRIntegrator()->DropAllChannels(false);
      m_fsrmode=0;
    }
    if (m_fsrmode==2) {
      p_pshandler->FSRIntegrator()->DropAllChannels();
      p_pshandler->FSRIntegrator()->Add(p_fsrinterface);
      p_fsrinterface->SetAlpha(1.0);
      p_fsrinterface->SetAlphaSave(1.0);
      m_fsrmode=1;
    }
  }
}

void Semihard_QCD::CreateISRChannels() 
{
  PHASIC::Multi_Channel *isr=p_pshandler->ISRIntegrator();
  isr->DropAllChannels();
  double mass=sqrt(p_isrhandler->SprimeMin());
  PHASIC::Single_Channel *channel=NULL;
  if (m_xsecs.size()>0 && 
      m_xsecs[0]->Regulator()!=NULL && 
      m_xsecs[0]->Regulator()->Type()==PHASIC::rf::qcd_trivial) {
    mass=2.0*m_xsecs[0]->Regulator()->Parameters()[0]; 
    channel = 
      new PHASIC::Flat_ISR_V(1.0," isr",p_pshandler->GetInfo());
  }
  else {
    channel = 
      new PHASIC::Simple_Pole_Uniform_V(1.0," isr",p_pshandler->GetInfo());
  }
  channel->SetAlpha(1.0);
  channel->SetAlphaSave(1.0);
  isr->Add(channel);
  isr->Reset();
}

void Semihard_QCD::InitIntegrators() 
{
  if (p_isrhandler->On()) {
    p_isrhandler->SetPartonMasses(p_flavours);
    for (unsigned int i=0;i<m_xsecs.size();i++) 
      m_xsecs[i]->SetISR(p_isrhandler);
    p_pshandler->ISRIntegrator()->DropAllChannels();
    if (p_pshandler->KMRZIntegrator()!=NULL) 
      p_pshandler->KMRZIntegrator()->DropAllChannels();
    if (p_pshandler->KMRKPIntegrator()!=NULL) 
      p_pshandler->KMRKPIntegrator()->DropAllChannels();
  }
  p_pshandler->CreateIntegrators();
  p_pshandler->ISRIntegrator()->Reset();
  if (p_isrhandler->KMROn()) {
    p_pshandler->KMRZIntegrator()->Reset();
    p_pshandler->KMRKPIntegrator()->Reset();
  }
  p_pshandler->FSRIntegrator()->Reset();
  m_channels=true;
}
      
