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
  Multi_Channel(name), p_psh(psh), m_keyid("BEAM"),
  p_beamspectra(p_psh->GetBeamSpectra())
{
  m_beammode = p_beamspectra?p_beamspectra->BeamMode():BEAM::beammode::unknown;
  for (size_t i=0;i<2;i++)
    m_beamtype[i] = (p_beamspectra?
		     p_beamspectra->GetBeam(i)->Type():
		     BEAM::beamspectrum::unknown);
  for (double yexp=-.999;yexp<=1.0;yexp+=.999) m_yexponents.insert(yexp);
}

bool Beam_Channels::Initialize()
{
  //msg_Out()<<METHOD<<" for mode = "<<m_beammode<<" "
  //	   <<"and "<<m_beamparams.size()<<" parameters.\n";
  return MakeChannels();
}

bool Beam_Channels::MakeChannels()
{
  if (m_beamparams.size()>0) return CreateChannels();
  switch (m_beammode) {
  case beammode::relic_density:
    m_beamparams.push_back(Channel_Info(channel_type::simple,0.5));
    m_beamparams.push_back(Channel_Info(channel_type::simple,2.0));
    CheckForStructuresFromME();
    break;
  case beammode::DM_annihilation:
    m_beamparams.push_back(Channel_Info(channel_type::simple,3));
    CheckForStructuresFromME();
    break;
  case beammode::collider:
    if (!DefineColliderChannels()) {
      msg_Error()<<"Error in "<<METHOD<<" for collider set-up:\n"
		 <<"   Don't know how to deal with combination of beamspectra: "
		 <<m_beamtype[0]<<" + "<<m_beamtype[1]<<".\n"
		 <<"   Will not initialize integration over spectra.\n";
    }
    break;
  case beammode::unknown:
  default:
    msg_Error()<<"Error in "<<METHOD<<":\n"
	       <<"   Unknown beam type.\n"
	       <<"   Will not initialize integration over spectra.\n";
    return false;
  }
  return CreateChannels();
}

bool Beam_Channels::DefineColliderChannels() {
  // default collider setup - no spectra
  if (m_beamtype[0]==beamspectrum::monochromatic &&
      m_beamtype[1]==beamspectrum::monochromatic) return true;
  // one or two laser backscattering spectra with monochromatic beams
  if ((m_beamtype[0]==beamspectrum::monochromatic &&
       (m_beamtype[1]==beamspectrum::laser_backscattering ||
	m_beamtype[1]==beamspectrum::simple_Compton)) ||
      ((m_beamtype[0]==beamspectrum::laser_backscattering ||
	m_beamtype[0]==beamspectrum::simple_Compton) &&
       m_beamtype[1]==beamspectrum::monochromatic) ||
      ((m_beamtype[0]==beamspectrum::laser_backscattering ||
	m_beamtype[0]==beamspectrum::simple_Compton) &&
       (m_beamtype[1]==beamspectrum::laser_backscattering ||
	m_beamtype[1]==beamspectrum::simple_Compton))) {
    m_beamparams.push_back(Channel_Info(channel_type::laserback,
					p_beamspectra->Peak(),1.));
    CheckForStructuresFromME();
    return true;
  }
  // one or two EPA spectra with monochromatic beams
  // currently our EPA is completely collinear, with real photons:
  // - todo: add proper EPA, with virtual photons and a physical deflection angle of
  //         the emitters.
  if ((m_beamtype[0]==beamspectrum::monochromatic &&
       m_beamtype[1]==beamspectrum::EPA) ||
      (m_beamtype[0]==beamspectrum::EPA &&
       m_beamtype[1]==beamspectrum::monochromatic) ||
      (m_beamtype[0]==beamspectrum::EPA &&
       m_beamtype[1]==beamspectrum::EPA)) {
    double exponent = (int(m_beamtype[0]==beamspectrum::EPA)+
		       int(m_beamtype[1]==beamspectrum::EPA))*0.75;
    m_beamparams.push_back(Channel_Info(channel_type::simple,exponent));
    CheckForStructuresFromME();
    return true;
  }
  if (m_beamtype[0]==beamspectrum::spectrum_reader ||
      m_beamtype[1]==beamspectrum::spectrum_reader) {
    msg_Error()<<"Warning in "<<METHOD<<":\n"
	       <<"   Beam spectra from spectrum reader - "
	       <<"will have to find a way to parse relevant information.\n"
	       <<"   Will pretend  a simple pole is good enough.\n";
    m_beamparams.push_back(Channel_Info(channel_type::simple,0.5));
    return true;
  }
}

void Beam_Channels::CheckForStructuresFromME() {
  if (!p_psh->Process()) {
    msg_Error()<<"Warning in "<<METHOD<<":\n"
	       <<"   Phase space handler has no process information.\n"
	       <<"   This looks like a potential bug, will exit.\n";
    THROW(fatal_error,"No process information in phase space handler.")
  }
  size_t nfsrchannels = p_psh->FSRIntegrator()->Number();
  std::vector<int>    types(nfsrchannels,0);
  std::vector<double> masses(nfsrchannels,0.0), widths(nfsrchannels,0.0);
  bool onshellresonance(false), fromFSR(false);  
  for (size_t i=0;i<nfsrchannels;i++) {
    p_psh->FSRIntegrator()->ISRInfo(i,types[i],masses[i],widths[i]);
    channel_type::code type = channel_type::code(abs(types[i]));
    switch (type) {
    case channel_type::simple:
    case channel_type::leadinglog:
    case channel_type::laserback:
      continue;
    case channel_type::threshold:
      if (ATOOLS::IsZero(masses[i])) continue;
      fromFSR = true;
      m_beamparams.push_back(Channel_Info(type,masses[i],2.));
      break;
    case channel_type::resonance:
      if (ATOOLS::IsZero(masses[i])) continue;
      if (types[i]==-1) {
	p_psh->SetOSMass(masses[i]);
	onshellresonance = true;
      }
      fromFSR = true;
      m_beamparams.push_back(Channel_Info(type,masses[i],widths[i]));
      break;
    }
  }
  if (fromFSR) return;
  Flavour_Vector * resonances = p_psh->Process()->Process()->Resonances();
  if (resonances && !resonances->empty()) {
    for (size_t i=0;i<resonances->size();i++) {
      Flavour flav = (*resonances)[i];
      double mass = flav.Mass();
      if (ATOOLS::IsZero(mass)) continue;
      m_beamparams.push_back(Channel_Info(channel_type::resonance,
					  mass,flav.Width()));
    }
  }
}

bool Beam_Channels::CreateChannels()
{
  if (m_beamparams.size() < 1) return 0;
  int beams = int(p_beamspectra->ColliderMode());
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case channel_type::simple:
      AddSimplePole(i,beams);
      break;
    case channel_type::resonance:
      AddResonance(i,beams);
      break;
    case channel_type::threshold:
      AddThreshold(i,beams);
      break;
    case channel_type::laserback:
      AddLaserBackscattering(i,beams);
      break;
    case channel_type::leadinglog:
    case channel_type::unknown:
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   tried to construct channel for unknown type.\n"
		 <<"   Will ignore this channel and hope for the best.\n";
    }
  }
  //msg_Out()<<METHOD<<" created "<<channels.size()<<" channels:\n";
  //for (size_t i=0;i<channels.size();i++)
  //  msg_Out()<<"  "<<channels[i]->Name()<<" : "<<channels[i]->Alpha()<<"\n";
  //msg_Out()<<"----------------------------------------------\n";
  return 1;
}

void Beam_Channels::AddSimplePole(const size_t & chno,const int & beams) {
  if (m_beammode==beammode::relic_density) {
    Add(new Simple_Pole_RelicDensity(m_beamparams[chno].parameters[0],
				     m_keyid,p_psh->GetInfo()));
    return;
  }
  else if (m_beammode==beammode::DM_annihilation) {
    double mass1 = p_beamspectra->GetBeam(0)->Beam().Mass();
    double mass2 = p_beamspectra->GetBeam(1)->Beam().Mass();
    Add(new Simple_Pole_DM_Annihilation(m_beamparams[chno].parameters[0],
					mass1,mass2,m_keyid,p_psh->GetInfo()));
    return;
  }
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new Simple_Pole_Uniform(m_beamparams[chno].parameters[0],
				  m_keyid,p_psh->GetInfo(),beams));
      Add(new Simple_Pole_Central(m_beamparams[chno].parameters[0],
				  m_keyid,p_psh->GetInfo(),beams));
    }
    else if (beams==3) {
      Add(new Simple_Pole_Forward(m_beamparams[chno].parameters[0],(*yit),
				  m_keyid,p_psh->GetInfo()));
      Add(new Simple_Pole_Backward(m_beamparams[chno].parameters[0],(*yit),
				   m_keyid,p_psh->GetInfo()));
    }
  }
}


void Beam_Channels::AddResonance(const size_t & chno,const int & beams) {
  if (m_beammode==beammode::relic_density) {
    Add(new Resonance_RelicDensity(m_beamparams[chno].parameters[0],
				   m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo()));
    return;
  }
  /*
  else if (m_beammode==beammode::DM_annihilation) {
    double mass1 = p_beamspectra->GetBeam(0)->Beam().Mass();
    double mass2 = p_beamspectra->GetBeam(1)->Beam().Mass();
    Add(new Resonance_DM_Annihilation(m_beamparams[chno].parameters[0],
					mass1,mass2,m_keyid,p_psh->GetInfo()));
    return;
  }
  */
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new Resonance_Uniform(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],
				m_keyid,p_psh->GetInfo(),beams));
      Add(new Resonance_Central(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],
				m_keyid,p_psh->GetInfo(),beams));
    }
    else if (beams==3) {
      Add(new Resonance_Forward(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],
				(*yit),m_keyid,p_psh->GetInfo()));
      Add(new Resonance_Backward(m_beamparams[chno].parameters[0],
				 m_beamparams[chno].parameters[1],
				 (*yit),m_keyid,p_psh->GetInfo()));
    }
  }
}

void Beam_Channels::AddThreshold(const size_t & chno,const int & beams) {
  if (m_beammode==beammode::relic_density) return;
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new Threshold_Uniform(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo(),beams));
      Add(new Threshold_Central(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],m_keyid,p_psh->GetInfo(),beams));
    }
    else if (beams==3) {
      Add(new Threshold_Forward(m_beamparams[chno].parameters[0],
				m_beamparams[chno].parameters[1],
				(*yit),m_keyid,p_psh->GetInfo()));
      Add(new Threshold_Backward(m_beamparams[chno].parameters[0],
				 m_beamparams[chno].parameters[1],
				 (*yit),m_keyid,p_psh->GetInfo()));
    }
  }
}

void Beam_Channels::AddLaserBackscattering(const size_t & chno,const int & beams) {
  if (m_beammode==beammode::relic_density) return;
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new LBS_Compton_Peak_Uniform(m_beamparams[chno].parameters[1],
				       m_beamparams[chno].parameters[0],
				       m_keyid,p_psh->GetInfo(),beams));
      Add(new LBS_Compton_Peak_Central(m_beamparams[chno].parameters[1],
				       m_beamparams[chno].parameters[0],
				       m_keyid,p_psh->GetInfo(),beams));
    }
    else if (beams==3) {
      Add(new LBS_Compton_Peak_Forward(m_beamparams[chno].parameters[1],
				       m_beamparams[chno].parameters[0],
				       (*yit),m_keyid,p_psh->GetInfo()));
      Add(new LBS_Compton_Peak_Backward(m_beamparams[chno].parameters[1],
					m_beamparams[chno].parameters[0],
				        (*yit),m_keyid,p_psh->GetInfo()));
    }
  }
}
