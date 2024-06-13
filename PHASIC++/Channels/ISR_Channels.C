#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/ISR_Channel_Base.H"
#include "PHASIC++/Channels/Simple_Pole_Channels.H"
#include "PHASIC++/Channels/Resonance_Channels.H"
#include "PHASIC++/Channels/Threshold_Channels.H"
#include "PHASIC++/Channels/Leading_Log_Channels.H"
#include "PHASIC++/Channels/LBS_Compton_Peak_Channels.H"

using namespace PHASIC;
using namespace ATOOLS;

ISR_Channels::ISR_Channels(Phase_Space_Handler *const psh,
			     const std::string &name) :
  Multi_Channel(name), p_psh(psh), m_keyid("ISR"),
  p_isrhandler(p_psh->GetISRHandler()),
  p_yfshandler(p_psh->GetYFSHandler()),
  m_isrmode(p_isrhandler->Mode())
{
  if(p_yfshandler->Mode()==YFS::yfsmode::isr){
    m_isrmode = PDF::isrmode::lepton_lepton;
    for (size_t i=0;i<2;i++) m_isrtype[i] = PDF::isrtype::yfs;
    for (double spexp=0.5;spexp<=1.5;spexp+=0.5) m_spexponents.insert(spexp);

  }
  else{
    for (size_t i=0;i<2;i++)
      m_isrtype[i] = (p_isrhandler?p_isrhandler->Type(i):PDF::isrtype::intact);
    for (double yexp=-.999;yexp<=1.0;yexp+=.999) m_yexponents.insert(yexp);
    for (double spexp=0.5;spexp<=1.5;spexp+=0.5) m_spexponents.insert(spexp);
  }
}

bool ISR_Channels::MakeChannels()
{
  if (m_isrparams.size()>0) return CreateChannels();
  switch (m_isrmode) {
  case PDF::isrmode::none:
    return true;
  case PDF::isrmode::hadron_hadron:
    for (std::set<double>::iterator spit=m_spexponents.begin();
	 spit!=m_spexponents.end();spit++) {
      for (std::set<double>::iterator yit=m_yexponents.begin();
	   yit!=m_yexponents.end();yit++) {
	m_isrparams.push_back(Channel_Info(channel_type::simple,(*spit),(*yit)));
      }
    }
    break;
  case PDF::isrmode::lepton_hadron:
  case PDF::isrmode::hadron_lepton: 
    m_isrparams.push_back(Channel_Info(channel_type::simple,1.,1.));
    break;
  case PDF::isrmode::lepton_lepton:
    m_isrparams.push_back(Channel_Info(channel_type::leadinglog,
				       p_psh->Process()->ISR()->Exponent(1),
				       1.00000001,1.));
    break;
  case PDF::isrmode::unknown:
  default:
    msg_Error()<<"Error in "<<METHOD<<": unknown isr mode.\n"
	       <<"   Continue without channels and hope for the best.\n";
    return true;
  }
  CheckForStructuresFromME();
  return CreateChannels();
}

void ISR_Channels::CheckForStructuresFromME() {
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
      for (std::set<double>::iterator yit=m_yexponents.begin();
	   yit!=m_yexponents.end();yit++) {
	m_isrparams.push_back(Channel_Info(type,masses[i],2.,(*yit)));
      }
      if(p_yfshandler->Mode()!=YFS::yfsmode::off){
        m_isrparams.push_back(Channel_Info(type,masses[i],2.));
      }
      break;
    case channel_type::resonance:
      if (ATOOLS::IsZero(masses[i])) continue;
      if (ATOOLS::IsZero(widths[i])) continue;
      if (types[i]==-1) {
	p_psh->SetOSMass(masses[i]);
	onshellresonance = true;
      }
      fromFSR = true;
      for (std::set<double>::iterator yit=m_yexponents.begin();
	   yit!=m_yexponents.end();yit++) {
	m_isrparams.push_back(Channel_Info(type,masses[i],widths[i],(*yit)));
      }
      if(p_yfshandler->Mode()!=YFS::yfsmode::off){
        m_isrparams.push_back(Channel_Info(type,masses[i],widths[i]));
      }
      break;
    case channel_type::simple:
      ifp_yfshandler->Mode()!=YFS::yfsmode::off{
        m_isrparams.push_back(Channel_Info(type,1.,1));
      }
      break;
    case channel_type::leadinglog:
    case channel_type::laserback:
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
      for (std::set<double>::iterator yit=m_yexponents.begin();
	   yit!=m_yexponents.end();yit++) {
	m_isrparams.push_back(Channel_Info(channel_type::resonance,
					   mass,flav.Width(),(*yit)));
      }
    }
  }
}

bool ISR_Channels::CreateChannels()
{
  size_t collmode =
    1*size_t(m_isrtype[0]!=PDF::isrtype::intact &&
	     m_isrtype[0]!=PDF::isrtype::unknown) +
    2*size_t(m_isrtype[1]!=PDF::isrtype::intact &&
	     m_isrtype[1]!=PDF::isrtype::unknown);
  if(m_isrtype[0]==PDF::isrtype::yfs) collmode = 4;
  if(p_yfshandler->Mode()==YFS::yfsmode::isr) collmode = 4;
  if (m_isrparams.size() < 1 || collmode==0) return 0;
  for (size_t i=0;i<m_isrparams.size();i++) {
    switch (m_isrparams[i].type) {
    case channel_type::simple:
      AddSimplePole(i,collmode);
      break;
    case channel_type::resonance:
      AddResonance(i,collmode);
      break;
    case channel_type::threshold:
      AddResonance(i,collmode);
      break;
    case channel_type::leadinglog:
      AddLeadingLog(i,collmode);
      break;
    case channel_type::exponential:
    case channel_type::laserback:
    case channel_type::unknown:
    default:
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   tried to construct channel for unknown type.\n"
		 <<"   Will ignore this channel and hope for the best.\n";
    }
  }
  return true;
}

void ISR_Channels::AddSimplePole(const size_t & chno,const size_t & mode) {
  double spexp = m_isrparams[chno].parameters[0];
  double yexp  = m_isrparams[chno].parameters.size()>1?m_isrparams[chno].parameters[1]:0.;
  if (mode==3 && (m_isrmode==PDF::isrmode::hadron_hadron ||
		  m_isrmode==PDF::isrmode::lepton_lepton)) {
    if (dabs(yexp)<1.e-3) {
      Add(new Simple_Pole_Uniform(spexp,m_keyid,p_psh->GetInfo(),mode));
      Add(new Simple_Pole_Central(spexp,m_keyid,p_psh->GetInfo(),mode));
    }
    else {
      Add(new Simple_Pole_Forward(spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
      Add(new Simple_Pole_Backward(spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
    }
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::lepton_hadron) {
    Add(new Simple_Pole_Uniform(spexp,m_keyid,p_psh->GetInfo(),mode));
    Add(new Simple_Pole_Forward(spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::hadron_lepton) {
    Add(new Simple_Pole_Uniform(spexp,m_keyid,p_psh->GetInfo(),mode));
    Add(new Simple_Pole_Backward(spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==4){
    // YFS Channels
    Add(new Simple_Pole_YFS(spexp, m_keyid, p_psh->GetInfo()));
  }
  else {
    Add(new Simple_Pole_Central(spexp,m_keyid,p_psh->GetInfo(),mode));
  }
}

void ISR_Channels::AddResonance(const size_t & chno,const size_t & mode) {
  double mass  = m_isrparams[chno].parameters[0];
  double width = m_isrparams[chno].parameters[1];
  double yexp  = m_isrparams[chno].parameters.size()>2?m_isrparams[chno].parameters[2]:0.;
  if (mode==3 && (m_isrmode==PDF::isrmode::hadron_hadron ||
		  m_isrmode==PDF::isrmode::lepton_lepton)) {
    if (dabs(yexp)<1.e-3) {
      Add(new Resonance_Uniform(mass,width,m_keyid,p_psh->GetInfo(),mode));
      Add(new Resonance_Central(mass,width,m_keyid,p_psh->GetInfo(),mode));
    }
    else {
      if (p_psh->GetBeamSpectra()->GetBeam(1)->Bunch() != Flavour(kf_photon) || !IsZero(yexp + 0.999))
        Add(new Resonance_Forward(mass,width,yexp,m_keyid,p_psh->GetInfo(),mode));
      if (p_psh->GetBeamSpectra()->GetBeam(0)->Bunch() != Flavour(kf_photon) || !IsZero(yexp - 0.999))
        Add(new Resonance_Backward(mass,width,yexp,m_keyid,p_psh->GetInfo(),mode));
    }
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::lepton_hadron) {
    Add(new Resonance_Uniform(mass,width,m_keyid,p_psh->GetInfo(),mode));
    Add(new Resonance_Forward(mass,width,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::hadron_lepton) {
    Add(new Resonance_Uniform(mass,width,m_keyid,p_psh->GetInfo(),mode));
    Add(new Resonance_Backward(mass,width,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==4 && m_isrmode==PDF::isrmode::lepton_lepton) {
    Add(new Resonance_YFS(mass,width,m_keyid,p_psh->GetInfo()));
  }
  else {
    Add(new Resonance_Central(mass,width,m_keyid,p_psh->GetInfo(),mode));
  }
}

void ISR_Channels::AddThreshold(const size_t & chno,const size_t & mode) {
  double mass  = m_isrparams[chno].parameters[0];
  double spexp = m_isrparams[chno].parameters[1];
  double yexp  = m_isrparams[chno].parameters.size()>2?m_isrparams[chno].parameters[2]:0.;
  if (mode==3 && (m_isrmode==PDF::isrmode::hadron_hadron ||
	       m_isrmode==PDF::isrmode::lepton_lepton)) {
    if (yexp==0.0) {
      Add(new Threshold_Uniform(mass,spexp,m_keyid,p_psh->GetInfo(),mode));
      Add(new Threshold_Central(mass,spexp,m_keyid,p_psh->GetInfo(),mode));
    }
    else {
      Add(new Threshold_Forward(mass,spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
      Add(new Threshold_Backward(mass,spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
    }
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::lepton_hadron) {
    Add(new Threshold_Uniform(mass,spexp,m_keyid,p_psh->GetInfo(),mode));
    Add(new Threshold_Forward(mass,spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::hadron_lepton) {
    Add(new Threshold_Uniform(mass,spexp,m_keyid,p_psh->GetInfo(),mode));
    Add(new Threshold_Backward(mass,spexp,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==4 && m_isrmode==PDF::isrmode::lepton_lepton) {
    Add(new Threshold_YFS(mass,spexp,m_keyid,p_psh->GetInfo()));
  }
  else {
    Add(new Threshold_Central(mass,spexp,m_keyid,p_psh->GetInfo(),mode));
  }
}

void ISR_Channels::AddLeadingLog(const size_t & chno,const size_t & mode) {
  double beta   = m_isrparams[chno].parameters[0];
  double factor = m_isrparams[chno].parameters[1];
  double yexp   = m_isrparams[chno].parameters.size()>2?m_isrparams[chno].parameters[2]:0.;
  if (mode==3 && m_isrmode==PDF::isrmode::lepton_lepton) {
    if (yexp==0.0) {
      Add(new Leading_Log_Uniform(beta,factor,m_keyid,p_psh->GetInfo(),mode));
      Add(new Leading_Log_Central(beta,factor,m_keyid,p_psh->GetInfo(),mode));
    }
    else {
      Add(new Leading_Log_Forward(beta,factor,yexp,m_keyid,p_psh->GetInfo(),mode));
      Add(new Leading_Log_Backward(beta,factor,yexp,m_keyid,p_psh->GetInfo(),mode));
    }
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::lepton_hadron) {
    Add(new Leading_Log_Uniform(beta,factor,m_keyid,p_psh->GetInfo(),mode));
    Add(new Leading_Log_Forward(beta,factor,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==3 && m_isrmode==PDF::isrmode::hadron_lepton) {
    Add(new Leading_Log_Uniform(beta,factor,m_keyid,p_psh->GetInfo(),mode));
    Add(new Leading_Log_Backward(beta,factor,yexp,m_keyid,p_psh->GetInfo(),mode));
  }
  else if (mode==4 && m_isrmode==PDF::isrmode::lepton_lepton) {
    Add(new Leading_Log_YFS(beta,factor,m_keyid,p_psh->GetInfo()));
  }
  else {
    Add(new Leading_Log_Uniform(beta,factor,m_keyid,p_psh->GetInfo(),mode));
    Add(new Leading_Log_Central(beta,factor,m_keyid,p_psh->GetInfo(),mode));
  }
}

