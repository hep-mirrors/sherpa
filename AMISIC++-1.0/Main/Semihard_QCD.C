#include "Semihard_QCD.H"

#include "FSR_Channel.H"
#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"

#define INIFL ATOOLS::Flavour(ATOOLS::kf::jet)

using namespace AMISIC;

Semihard_QCD::Semihard_QCD(BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr,
			   ATOOLS::Selector_Data *const seldata,
			   const ATOOLS::Flavour *flavours,
			   const int scalescheme,
			   const int kfactorscheme,
			   const double scalefactor):
  XS_Group(2,2,flavours,scalescheme,kfactorscheme,scalefactor,
	   beam,isr,seldata)
{
  SetFSRInterface(NULL);
  SetFSRMode(0);
}

Semihard_QCD::~Semihard_QCD()
{
}

void Semihard_QCD::CreateFSRChannels() 
{
  if ((m_fsrmode==0)||(p_fsrinterface==NULL)) {
    p_pshandler->FSRIntegrator()->DropAllChannels();
    p_pshandler->FSRIntegrator()->
      Add(new PHASIC::S1Channel(2,2,p_flavours,
				ATOOLS::Flavour(ATOOLS::kf::gluon)));
    p_pshandler->FSRIntegrator()->Add(new PHASIC::T1Channel(2,2,p_flavours));
    p_pshandler->FSRIntegrator()->Add(new PHASIC::U1Channel(2,2,p_flavours));
    m_fsrmode=1;
  }
  else {
    if (m_fsrmode==2) {
      p_pshandler->FSRIntegrator()->DropAllChannels();
      p_pshandler->FSRIntegrator()->Add(p_fsrinterface);
      m_fsrmode=1;
    }
  }
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
}
      
