#include "PHASIC++/Main/Channel_Creator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Beam_Channels.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/CXXFLAGS.H"

using namespace PHASIC;
using namespace ATOOLS;

Channel_Creator::Channel_Creator(Phase_Space_Handler * psh) :
  p_psh(psh) {}

Channel_Creator::~Channel_Creator() {}

bool Channel_Creator::operator()() 
{
  msg_Tracking()<<"Initializing channels for phase space integration (\n\t";
  Settings& s = Settings::GetMainSettings();
  if (!CreateFSRIntegrator())  THROW(fatal_error,"Could not create FSR channels");
  if (!CreateBeamIntegrator()) THROW(fatal_error,"Could not create beam channels");
  if (!CreateISRIntegrator())  THROW(fatal_error,"Could not create ISR channels");
  msg_Tracking()<<")\n";
  return true;
}

bool Channel_Creator::CreateBeamIntegrator() {
  if (p_psh->Process()->NIn()!=2) return true;
  BEAM::Beam_Spectra_Handler * beamhandler = p_psh->GetBeamSpectra();
  if (!beamhandler || !beamhandler->On()) return true;
  Beam_Channels * beamchannels =
    new Beam_Channels(p_psh,"beam_"+p_psh->Process()->Process()->Name());
  p_psh->SetBeamIntegrator(beamchannels);
  return beamchannels->Initialize();
}

bool Channel_Creator::CreateISRIntegrator() {
  if (p_psh->Process()->NIn()!=2) return true;
  PDF::ISR_Handler * isrhandler = p_psh->GetISRHandler();
  if (!isrhandler || isrhandler->On()==0) return true;
  ISR_Channels * isrchannels =
    new ISR_Channels(p_psh,"isr_"+p_psh->Process()->Process()->Name());
  p_psh->SetISRIntegrator(isrchannels);
  return isrchannels->Initialize();
}

bool Channel_Creator::CreateFSRIntegrator() {
  FSR_Channels * fsrchannels =
    new FSR_Channels(p_psh,"fsr_"+p_psh->Process()->Process()->Name());
  p_psh->SetFSRIntegrator(fsrchannels);
  return fsrchannels->Initialize();
}

