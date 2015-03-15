#include "PHASIC++/Channels/ISR_Channels.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Channels/ISR_Vegas.H"

using namespace PHASIC;
using namespace ATOOLS;

bool ISR_Channels::MakeChannels()
{
  if (m_isrparams.size()>0) return CreateChannels();
  bool onshellresonance=false;
  Channel_Info ci;
  int    type;
  double mass,width;
  std::set<double> ths;
  ths.insert(sqrt(p_psh->Cuts()->Smin()));
  Multi_Channel *fsr(p_psh->FSRIntegrator());
  std::vector<int> ts(fsr->Number(),0);
  std::vector<double> ms(fsr->Number(),0.0), ws(fsr->Number(),0.0);
  for (size_t i(0);i<fsr->Number();++i) fsr->ISRInfo(i,ts[i],ms[i],ws[i]);
  fsr->ISRInfo(ts,ms,ws);
  for (size_t i=0;i<ts.size();i++) {
    type=abs(ts[i]); 
    if (ts[i]==-1) {
      p_psh->SetOSMass(ms[i]);
      onshellresonance=true;
    }
    mass=ms[i];
    width=ws[i];
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    if (type==2) {
      ths.insert(mass);
      continue;
    }
    for (double yexp=-.999;yexp<=1.0;yexp+=.999) {
    ci.type = type;
    (ci.parameters).push_back(mass);
    if (type==1) (ci.parameters).push_back(width);
    if (type==2) (ci.parameters).push_back(2.);
    if (p_psh->Flavs()[0].IsLepton() || 
	p_psh->Flavs()[1].IsLepton()) (ci.parameters).push_back(yexp);
    else (ci.parameters).push_back(yexp);
    bool add=true;
    for (size_t j=0;j<m_isrparams.size();j++) if (m_isrparams[j]==ci) add=false; 
    if (add) m_isrparams.push_back(ci);
    ci.parameters.clear();
    }
  }
  if (ths.size()) {
    for (std::set<double>::const_iterator thit(ths.begin());thit!=ths.end();++thit)
    for (double yexp=-.999;yexp<=1.0;yexp+=.999) {
    ci.type = 2;
    (ci.parameters).push_back(*thit);
    (ci.parameters).push_back(2.);
    if (p_psh->Flavs()[0].IsLepton() || 
	p_psh->Flavs()[1].IsLepton()) (ci.parameters).push_back(yexp);
    else (ci.parameters).push_back(yexp);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    }
  }
  if (onshellresonance) return CreateChannels();
  
  if (p_psh->Flavs()[0].IsLepton() || p_psh->Flavs()[1].IsLepton()) {
    if ((p_psh->Flavs()[0].IsLepton() && p_psh->Flavs()[1].Strong()) ||
	(p_psh->Flavs()[1].IsLepton() && p_psh->Flavs()[0].Strong())) {
      //The DIS case
      ci.type = 0; 
      (ci.parameters).push_back(1.);
      (ci.parameters).push_back(1.);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
    }
    else {
      // leptons : 1/s'^2 and 1/(s-s')^beta, sharp FW-BW peak
      //     ci.type = 0;
      //     (ci.parameters).push_back(.5);
      //     (ci.parameters).push_back(1.);
      //     m_isrparams.push_back(ci);
      //     ci.parameters.clear();
      //     ci.type = 0;
      //     (ci.parameters).push_back(2.);
      //     (ci.parameters).push_back(1.);
      //     m_isrparams.push_back(ci);
      //     ci.parameters.clear();
      ci.type = 3;
      (ci.parameters).push_back
	(p_psh->Process()->ISR()->Exponent(1));
      (ci.parameters).push_back(1.00000001);
      (ci.parameters).push_back(1.);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  else {
    // default : 1/s'
    for (double sexp=0.5;sexp<=1.5;sexp+=0.5) {
      for (double yexp=-.999;yexp<=1.0;yexp+=.999) {
	ci.type = 0;
	(ci.parameters).push_back(sexp);
	(ci.parameters).push_back(yexp);
	m_isrparams.push_back(ci);
	ci.parameters.clear();
      }
    }
  }
  return CreateChannels();
}

bool ISR_Channels::CreateChannels()
{
  if (m_isrparams.size() < 1) return 0;
  int isr = p_psh->Process()->ISR()->On();
  int ep = (p_psh->Flavs()[0].IsLepton() && p_psh->Flavs()[1].Strong()) ||
    (p_psh->Flavs()[1].IsLepton() && p_psh->Flavs()[0].Strong());
  for (size_t i=0;i<m_isrparams.size();i++) {
    switch (m_isrparams[i].type) {
    case 0:
      if (isr==3 && !ep) {
	if (m_isrparams[i].parameters[1]==0.0) {
 	Add(new Simple_Pole_Uniform_V
	    (m_isrparams[i].parameters[0]," isr",p_psh->GetInfo()));
	Add(new Simple_Pole_Central_V
	    (m_isrparams[i].parameters[0]," isr",p_psh->GetInfo(),isr));
	}
	else {
	Add(new Simple_Pole_Forward_V
	    (m_isrparams[i].parameters[0],
	     m_isrparams[i].parameters[1]," isr",p_psh->GetInfo()));
	Add(new Simple_Pole_Backward_V
	    (m_isrparams[i].parameters[0],
	     m_isrparams[i].parameters[1]," isr",p_psh->GetInfo()));
	}
      }
      else if (isr==3 && ep) {
	Add(new Simple_Pole_Forward_V
	    (m_isrparams[i].parameters[0],
	     m_isrparams[i].parameters[1]," isr",p_psh->GetInfo()));
      }
      else {
	Add(new Simple_Pole_Central_V
	    (m_isrparams[i].parameters[0]," isr",p_psh->GetInfo(),isr));
      }
      break;
    case 1:
      if (isr==3 && !ep) {
	if (m_isrparams[i].parameters[2]==0.0) {
	Add(new Resonance_Uniform_V
	    (m_isrparams[i].parameters[0],
	     m_isrparams[i].parameters[1]," isr",p_psh->GetInfo()));
	}
	else {
	Add(new Resonance_Forward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
	Add(new Resonance_Backward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
	}
      }
      else if (isr==3 && ep) {
	Add(new Resonance_Forward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
      }
      else {
	Add(new Resonance_Central_V
	    (m_isrparams[i].parameters[0],
	     m_isrparams[i].parameters[1]," isr",p_psh->GetInfo(),isr));
      }
      break;
    case 2:
      if (m_isrparams[i].parameters[2]==0.0) {
      Add(new Threshold_Central_V
	  (m_isrparams[i].parameters[0],
	   m_isrparams[i].parameters[1]," isr",p_psh->GetInfo(),isr));
      }
      else {
      if (isr==3) {
	Add(new Threshold_Forward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
	Add(new Threshold_Backward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
      }
      }
      break;
    case 3:
      Add(new Leading_Log_Central_V
	  (m_isrparams[i].parameters[0],
	   m_isrparams[i].parameters[1]," isr",p_psh->GetInfo(),isr));
      if (isr==3) {
	Add(new Leading_Log_Forward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
	Add(new Leading_Log_Backward_V
	    (m_isrparams[i].parameters[0],m_isrparams[i].parameters[1],
	     m_isrparams[i].parameters[2]," isr",p_psh->GetInfo()));
      }
      break;
    }
  }
  return 1;
}

bool ISR_Channels::Initialize()
{
  return MakeChannels();
}
