#include "PHASIC++/Channels/Beam_Channels.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/ISR_Channel_Base.H"
#include "PHASIC++/Channels/Simple_Pole_Channels.H"
#include "PHASIC++/Channels/Resonance_Channels.H"
#include "PHASIC++/Channels/Threshold_Channels.H"
#include "PHASIC++/Channels/Leading_Log_Channels.H"
#include "PHASIC++/Channels/LBS_Compton_Peak_Channels.H"

using namespace PHASIC;
using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Beam_Channels::Beam_Channels(Phase_Space_Handler *const psh,
		     const std::string &name) :
  Multi_Channel(name), p_psh(psh), m_keyid("BEAM") {
  BEAM::Beam_Spectra_Handler * beamspectra = p_psh->GetBeamSpectra();
  m_beammode = beamspectra?beamspectra->Mode():BEAM::beammode::unknown;
}

bool Beam_Channels::Initialize()
{
  msg_Out()<<METHOD<<" for mode = "<<m_beammode<<" "
  	   <<"and "<<m_beamparams.size()<<" parameters.\n";
  return MakeChannels();
}

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
  else if ((p_psh->Flavs()[0].IsPhoton()) || (p_psh->Flavs()[1].IsPhoton())) {
    // Laser Backscattering spectrum
    ci.type = 3;
    (ci.parameters).push_back(p_psh->Process()->Beam()->Peak());
    (ci.parameters).push_back(p_psh->Process()->Beam()->Exponent(1));
    (ci.parameters).push_back(0.7);
    m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  else {
    if (m_beammode!=beammode::relic_density) {
      ci.type = 0;
      (ci.parameters).push_back(.5);
      (ci.parameters).push_back(0.99);
      m_beamparams.push_back(ci);
      ci.parameters.clear();
    }
    if (m_beammode==beammode::relic_density) {
      ci.type = 0;
      (ci.parameters).push_back(.5);
      m_beamparams.push_back(ci);
      ci.parameters.clear();
      ci.type = 0;
      (ci.parameters).push_back(2.);
      m_beamparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  int    type;
  double mass,width;
  double thmin=0.,thmax=0.;
  msg_Out()<<"  trying to obtain resonacen information from FSR:\n";
  for (size_t i=0;i<p_psh->FSRIntegrator()->Number();i++) {
    type=0; 
    mass=width=0.;
    msg_Out()<<"   --> trying "<<i<<": ";
    if (p_psh->Process()) {
      p_psh->FSRIntegrator()->ISRInfo(i,type,mass,width);
      msg_Out()<<type<<", mass = "<<mass<<", width = "<<width<<"\n";
    }
    else msg_Out()<<"no process.\n";
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
    if ((p_psh->Flavs()[0].IsLepton()) || (p_psh->Flavs()[1].IsLepton()))
      (ci.parameters).push_back(1.);
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
    if ((p_psh->Flavs()[0].IsLepton()) || (p_psh->Flavs()[1].IsLepton()))
      (ci.parameters).push_back(1.);
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
  //msg_Out()<<METHOD<<" for "<<m_beamparams.size()<<" parameters, "
  //	   <<"beam = "<<p_psh->Process()->Beam()->On()<<" and info = "<<p_psh->GetInfo()<<"\n";
  if (m_beamparams.size() < 1) return 0;
  int beam = p_psh->Process()->Beam()->On();
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case 0:
      AddSimplePole(i,beam);
      break;
    case 1:
      AddResonance(i,beam);
      break;
    case 2:
      AddThreshold(i,beam);
      break;
    case 3:
      AddLaserBackscattering(i,beam);
      break;
    }
  }
  //msg_Out()<<METHOD<<" created "<<channels.size()<<" channels:\n";
  //for (size_t i=0;i<channels.size();i++) 
  //  msg_Out()<<"  "<<channels[i]->Name()<<" : "<<channels[i]->Alpha()<<"\n";
  //msg_Out()<<"----------------------------------------------\n";
  return 1;
}

void Beam_Channels::AddSimplePole(const size_t & chno,const int & beam) {
  if (m_beammode==beammode::relic_density) {
    Add(new Simple_Pole_RelicDensity(m_beamparams[chno].parameters[0],
				     m_keyid,p_psh->GetInfo()));
    return;
  }
  Add(new Simple_Pole_Central(m_beamparams[chno].parameters[0],
			      m_keyid,p_psh->GetInfo(),beam));
  if (beam!=3) return;
  Add(new Simple_Pole_Forward(m_beamparams[chno].parameters[0],
			      m_beamparams[chno].parameters[1],
			      m_keyid,p_psh->GetInfo()));
  Add(new Simple_Pole_Backward(m_beamparams[chno].parameters[0],
			       m_beamparams[chno].parameters[1],
			       m_keyid,p_psh->GetInfo()));
}


void Beam_Channels::AddResonance(const size_t & chno,const int & beam) {
  if (m_beammode==beammode::relic_density) {
    Add(new Resonance_RelicDensity(m_beamparams[chno].parameters[0],
				   m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo()));
    return;
  }
  Add(new Resonance_Central(m_beamparams[chno].parameters[0],
			    m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo(),beam));
  if (beam!=3) return;
  Add(new Resonance_Uniform(m_beamparams[chno].parameters[0],
			    m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo()));
  Add(new Resonance_Forward(m_beamparams[chno].parameters[0],
			    m_beamparams[chno].parameters[1],
			    m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
  Add(new Resonance_Backward(m_beamparams[chno].parameters[0],
			     m_beamparams[chno].parameters[1],
			     m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
}
  
void Beam_Channels::AddThreshold(const size_t & chno,const int & beam) {
  if (m_beammode==beammode::relic_density) return;
  Add(new Threshold_Central(m_beamparams[chno].parameters[0],
			    m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo(),beam));
  if (beam!=3) return;
  Add(new Threshold_Forward(m_beamparams[chno].parameters[0],
			    m_beamparams[chno].parameters[1],
			    m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
  Add(new Threshold_Backward(m_beamparams[chno].parameters[0],
			     m_beamparams[chno].parameters[1],
			     m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
}

void Beam_Channels::AddLaserBackscattering(const size_t & chno,const int & beam) {
  if (m_beammode==beammode::relic_density) return;
  if (!p_psh->Flavs()[0].IsPhoton() && !p_psh->Flavs()[1].IsPhoton()) return;
  Add(new LBS_Compton_Peak_Central(m_beamparams[chno].parameters[1],
				   m_beamparams[chno].parameters[0],m_keyid,p_psh->GetInfo(),beam));
  if (beam!=3) return;
  Add(new LBS_Compton_Peak_Forward(m_beamparams[chno].parameters[1],
				   m_beamparams[chno].parameters[0],
				   m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
  Add(new LBS_Compton_Peak_Backward(m_beamparams[chno].parameters[1],
				    m_beamparams[chno].parameters[0],
				    m_beamparams[chno].parameters[2],m_keyid,p_psh->GetInfo()));
}

