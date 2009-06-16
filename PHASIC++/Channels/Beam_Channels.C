#include "PHASIC++/Channels/Beam_Channels.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/ISR_Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

bool Beam_Channels::MakeChannels()
{
  if (m_beamparams.size()>0) return CreateChannels();
  Channel_Info ci;
  // default : Beamstrahlung
  if ((p_psh->Flavs()[0].IsLepton()) && (p_psh->Flavs()[1].IsLepton())) {
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(p_psh->Process()->Beam()->Exponent(1));
    m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  else {
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(0.99);
    m_beamparams.push_back(ci);
    ci.parameters.clear();
    // Laser Backscattering spectrum
    ci.type = 3;
    (ci.parameters).push_back(p_psh->Process()->Beam()->Peak());
    (ci.parameters).push_back(p_psh->Process()->Beam()->Exponent(1));
    (ci.parameters).push_back(0.7);
    m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  int    type;
  double mass,width;
  double thmin=0.,thmax=0.;
  for (size_t i=0;i<p_psh->FSRIntegrator()->Number();i++) {
    type=0; 
    mass=width=0.;
    if (p_psh->Process()) p_psh->FSRIntegrator()->ISRInfo(i,type,mass,width);
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    if (type==2) {
      if (thmax==0.) { thmax=mass; thmin=mass; }
      thmin = ATOOLS::Min(thmin,mass);
      thmax = ATOOLS::Max(thmax,mass);
      continue;
    }
    ci.type = type;
    (ci.parameters).push_back(mass);
    if (type==1) (ci.parameters).push_back(width);
    if (type==2) (ci.parameters).push_back(1.5);
    if ((p_psh->Flavs()[0].IsLepton()) || (p_psh->Flavs()[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    bool add=true;
    for (size_t j=0;j<m_beamparams.size();j++) if (m_beamparams[j]==ci) add=false;
    if (add) m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  if (thmax>0.) {
    ci.type = 2;
    (ci.parameters).push_back(thmax);
    (ci.parameters).push_back(1.5);
    if ((p_psh->Flavs()[0].IsLepton()) || (p_psh->Flavs()[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    m_beamparams.push_back(ci);
    if (thmin<thmax) {
      (ci.parameters)[0]=thmin;
      m_beamparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  return CreateChannels();
}

bool Beam_Channels::CreateChannels()
{
  if (m_beamparams.size() < 1) return 0;
  int beam = p_psh->Process()->Beam()->On();
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case 0:
      Add(new Simple_Pole_Central_V(m_beamparams[i].parameters[0]," beam",p_psh->GetInfo(),beam));
      if (beam==3) {
	Add(new Simple_Pole_Forward_V(m_beamparams[i].parameters[0],
				      m_beamparams[i].parameters[1]," beam",p_psh->GetInfo()));
	Add(new Simple_Pole_Backward_V(m_beamparams[i].parameters[0],
				       m_beamparams[i].parameters[1]," beam",p_psh->GetInfo()));
      }
      break;
    case 1:
      Add(new Resonance_Central_V(m_beamparams[i].parameters[0],
				  m_beamparams[i].parameters[1]," beam",p_psh->GetInfo(),beam));
      if (beam==3) {
	Add(new Resonance_Uniform_V(m_beamparams[i].parameters[0],
				    m_beamparams[i].parameters[1]," beam",p_psh->GetInfo()));
	Add(new Resonance_Forward_V(m_beamparams[i].parameters[0],
				    m_beamparams[i].parameters[1],
				    m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
	Add(new Resonance_Backward_V(m_beamparams[i].parameters[0],
				     m_beamparams[i].parameters[1],
				     m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
      }
      break;
    case 2:
      Add(new Threshold_Central_V(m_beamparams[i].parameters[0],
				  m_beamparams[i].parameters[1]," beam",p_psh->GetInfo(),beam));
      if (beam==3) {
	Add(new Threshold_Forward_V(m_beamparams[i].parameters[0],
				    m_beamparams[i].parameters[1],
				    m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
	Add(new Threshold_Backward_V(m_beamparams[i].parameters[0],
				     m_beamparams[i].parameters[1],
				     m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
      }
      break;
    case 3:
      if ((p_psh->Flavs()[0].IsPhoton()) || (p_psh->Flavs()[1].IsPhoton())) {
	Add(new LBS_Compton_Peak_Central_V(m_beamparams[i].parameters[1],
					   m_beamparams[i].parameters[0]," beam",p_psh->GetInfo(),beam));
	if (beam==3) {
	  Add(new LBS_Compton_Peak_Forward_V(m_beamparams[i].parameters[1],
					     m_beamparams[i].parameters[0],
					     m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
	  Add(new LBS_Compton_Peak_Backward_V(m_beamparams[i].parameters[1],
					      m_beamparams[i].parameters[0],
					      m_beamparams[i].parameters[2]," beam",p_psh->GetInfo()));
	}
      }
      break;
    }
  }
  return 1;
}

bool Beam_Channels::Initialize()
{
  return MakeChannels();
}
