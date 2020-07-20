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
#include "PHASIC++/Channels/Exponential_Channels.H"

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
  return MakeChannels();
}

bool Beam_Channels::MakeChannels()
{
  if (m_beamparams.size()>0) return CreateChannels();
  switch (m_beammode) {
  case beammode::relic_density:
    {
      m_beamparams.push_back(Channel_Info(channel_type::simple,10.));
      m_beamparams.push_back(Channel_Info(channel_type::simple,1.));
    }
    CheckForStructuresFromME();
    break;
  case beammode::DM_annihilation:
    {
      Settings& settings = Settings::GetMainSettings();
      double temperature = settings["DM_TEMPERATURE"].Get<double>();
      double sexp = 1./(2*pow(temperature,2));
      m_beamparams.push_back(Channel_Info(channel_type::simple,1.));
      m_beamparams.push_back(Channel_Info(channel_type::exponential,sexp));
    }
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
					1.,p_beamspectra->Peak()));
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
  }
  return true;
}

void Beam_Channels::CheckForStructuresFromME() {
  if (!p_psh->Process()) {
    msg_Error()<<"Warning in "<<METHOD<<":\n"
	       <<"   Phase space handler has no process information.\n"
	       <<"   This looks like a potential bug, will exit.\n";
    THROW(fatal_error,"No process information in phase space handler.")
  }
  std::set<double>    thresholds;
  if (p_psh->Flavs()[0].Strong() && p_psh->Flavs()[1].Strong()) {
    if (p_psh->Cuts()!=NULL) thresholds.insert(sqrt(p_psh->Cuts()->Smin()));
  }
  size_t nfsrchannels = p_psh->FSRIntegrator()->Number();
  std::vector<int>    types(nfsrchannels,0);
  std::vector<double> masses(nfsrchannels,0.0), widths(nfsrchannels,0.0);
  for (size_t i=0;i<nfsrchannels;i++) {
    p_psh->FSRIntegrator()->ISRInfo(i,types[i],masses[i],widths[i]);
  }
  p_psh->FSRIntegrator()->ISRInfo(types,masses,widths);
  bool onshellresonance(false), fromFSR(false);  
  for (size_t i=0;i<types.size();i++) {    
    channel_type::code type = channel_type::code(abs(types[i]));
    switch (type) {
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
    case channel_type::simple:
    case channel_type::leadinglog:
    case channel_type::laserback:
		case channel_type::exponential:
    default:
      break;
    }
  }
  if (fromFSR) return;
  Flavour_Vector resonances;
  msg_Out()<<METHOD<<" for "<<fromFSR<<": "<<resonances<<"\n";
  if (p_psh->Process()->Process()->FillResonances(resonances) &&
      !resonances.empty()) {
    for (size_t i=0;i<resonances.size();i++) {
      Flavour flav = resonances[i];
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
  size_t mode = 0;
  if (p_beamspectra->BeamMode()==BEAM::beammode::collider) {
    switch (p_beamspectra->ColliderMode()) {
    case BEAM::collidermode::spectral_1:
      mode = 1;
      break;
    case BEAM::collidermode::spectral_2:
      mode = 2;
      break;
    case BEAM::collidermode::both_spectral:
      mode = 3;
      break;
    default:
      break;
    }
  }
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case channel_type::simple:
      AddSimplePole(i,mode);
      break;
    case channel_type::resonance:
      AddResonance(i,mode);
      break;
    case channel_type::threshold:
      AddThreshold(i,mode);
      break;
    case channel_type::laserback:
      AddLaserBackscattering(i,mode);
      break;
    case channel_type::exponential:
			AddExponential(i,mode);
      break;
    case channel_type::leadinglog:
    case channel_type::unknown:
    default:
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   tried to construct channel for unknown type.\n"
		 <<"   Will ignore this channel and hope for the best.\n";
    }
  }
  return true;
}

void Beam_Channels::AddSimplePole(const size_t & chno,const size_t & mode) {
  double spex = m_beamparams[chno].parameters[0];
  if (m_beammode==beammode::relic_density) {
    Add(new Simple_Pole_RelicDensity(spex,m_keyid,p_psh->GetInfo()));
  }
  else if (m_beammode==beammode::DM_annihilation) {
    double mass1 = p_beamspectra->GetBeam(0)->Beam().Mass();
    double mass2 = p_beamspectra->GetBeam(1)->Beam().Mass();
    Add(new Simple_Pole_DM_Annihilation(spex,mass1,mass2,m_keyid,p_psh->GetInfo()));
  }
  else {
    for (set<double>::iterator yit=m_yexponents.begin();
	 yit!=m_yexponents.end();yit++) {
      if (dabs(*yit)<1.e-3) {
	Add(new Simple_Pole_Uniform(spex,m_keyid,p_psh->GetInfo(),mode));
	Add(new Simple_Pole_Central(spex,m_keyid,p_psh->GetInfo(),mode));
      }
      else if (mode==3) {
	Add(new Simple_Pole_Forward(spex,(*yit),m_keyid,p_psh->GetInfo(),mode));
	Add(new Simple_Pole_Backward(spex,(*yit),m_keyid,p_psh->GetInfo(),mode));
      }
    }
  }
}

void Beam_Channels::AddResonance(const size_t & chno,const size_t & mode) {
  double mass  = m_beamparams[chno].parameters[0];
  double width = m_beamparams[chno].parameters[1];
  if (m_beammode==beammode::relic_density) {
    Add(new Resonance_RelicDensity(mass,
				   width,m_keyid,p_psh->GetInfo()));
    return;
  }
  else if (m_beammode==beammode::DM_annihilation) {
    double mass1 = p_beamspectra->GetBeam(0)->Beam().Mass();
    double mass2 = p_beamspectra->GetBeam(1)->Beam().Mass();
    Add(new Resonance_DM_Annihilation(mass,width,mass1,mass2,m_keyid,p_psh->GetInfo()));
    return;
  }
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new Resonance_Uniform(mass,width,m_keyid,p_psh->GetInfo(),mode));
      Add(new Resonance_Central(mass,width,m_keyid,p_psh->GetInfo(),mode));
    }
    else if (mode==3) {
      Add(new Resonance_Forward(mass,width,(*yit),m_keyid,p_psh->GetInfo(),mode));
      Add(new Resonance_Backward(mass,width,(*yit),m_keyid,p_psh->GetInfo(),mode));
    }
  }
}

void Beam_Channels::AddThreshold(const size_t & chno,const size_t & mode) {
  double mass = m_beamparams[chno].parameters[0];
  double spex = m_beamparams[chno].parameters[1];
  if (m_beammode==beammode::relic_density) return;
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new Threshold_Uniform(mass,spex,m_keyid,p_psh->GetInfo(),mode));
      Add(new Threshold_Central(mass,spex,m_keyid,p_psh->GetInfo(),mode));
    }
    else if (mode==3) {
      Add(new Threshold_Forward(mass,spex,(*yit),m_keyid,p_psh->GetInfo(),mode));
      Add(new Threshold_Backward(mass,spex,(*yit),m_keyid,p_psh->GetInfo(),mode));
    }
  }
}

void Beam_Channels::AddLaserBackscattering(const size_t & chno,const size_t & mode) {
  double exponent = m_beamparams[chno].parameters[0];
  double pole     = m_beamparams[chno].parameters[1];
  if (m_beammode==beammode::relic_density) return;
  for (set<double>::iterator yit=m_yexponents.begin();
       yit!=m_yexponents.end();yit++) {
    if (dabs(*yit)<1.e-3) {
      Add(new LBS_Compton_Peak_Uniform(exponent,pole,m_keyid,p_psh->GetInfo(),mode));
      Add(new LBS_Compton_Peak_Central(exponent,pole,m_keyid,p_psh->GetInfo(),mode));
    }
    else if (mode==3) {
      Add(new LBS_Compton_Peak_Forward(exponent,pole,(*yit),m_keyid,p_psh->GetInfo(),mode));
      Add(new LBS_Compton_Peak_Backward(exponent,pole,(*yit),m_keyid,p_psh->GetInfo(),mode));
    }
  }
}

void Beam_Channels::AddExponential(const size_t & chno,const size_t & mode) {
  double spex = m_beamparams[chno].parameters[0];
	double mass1 = p_beamspectra->GetBeam(0)->Beam().Mass();
	double mass2 = p_beamspectra->GetBeam(1)->Beam().Mass();
  if (m_beammode==beammode::relic_density) {
		Add(new Exponential_RelicDensity(spex,mass1,mass2,m_keyid,p_psh->GetInfo()));
  }
  else if (m_beammode==beammode::DM_annihilation) {
		// todo: change this
    Add(new Exponential_DM_Annihilation(spex,mass1,mass2,m_keyid,p_psh->GetInfo()));
  }
  else return;
}
