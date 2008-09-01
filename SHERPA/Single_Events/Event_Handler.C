#include "Event_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "My_Limits.H"
#include "Signal_Processes.H"
#include "CXXFLAGS.H"
#ifdef USING__PYTHIA
#include "Lund_Interface.H"
#endif
#include <unistd.h>
#include <cassert>

#include "Random.H"

#ifdef PROFILE__Event_Handler
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;


Event_Handler::Event_Handler():
  m_lastparticlecounter(0),
  m_lastblobcounter(0),
  p_mehandler(NULL)
{
  p_phases  = new Phase_List;
}

Event_Handler::~Event_Handler() 
{
  Reset();
  m_blobs.Clear();
  EmptyEventPhases();
  
  if (p_phases)   { delete p_phases;   p_phases   = NULL; }
}

void Event_Handler::AddEventPhase(Event_Phase_Handler * phase) 
{
  eph::code type   = phase->Type();
  std::string name = phase->Name();
  if (name.find("Signal_Processes:")!=std::string::npos) 
    p_mehandler=static_cast<Signal_Processes*>(phase)->GetMEHandler();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) { 
    if ((type==(*pit)->Type()) && (name==(*pit)->Name())) {
      msg_Out()<<"WARNING in Event_Handler::AddEventPhase("<<type<<":"<<name<<") "
		  <<"already included."<<std::endl;
      return;
    }
  }
  msg_Tracking()<<"Event_Handler::AddEventPhase("<<type<<":"<<name<<")."<<std::endl;
  p_phases->push_back(phase);
}

Event_Phase_Handler * Event_Handler::GetEventPhase(const size_t i) {
  if (i<p_phases->size()) {
    size_t count=i;
    for (Phase_Iterator pit=p_phases->begin();pit<p_phases->end();pit++) {
      if (count==0) return (*pit);
      count--;
    }
  }
  msg_Error()<<"Error in Event_Handler::GetEventPhase("<<i<<")"<<std::endl
	     <<"   Out of bounds, only "<<p_phases->size()<<" event phases."<<std::endl
	     <<"   Will return NULL."<<std::endl;
  return NULL;
}

size_t Event_Handler::NumberOfEventPhases() { return p_phases->size(); }

void Event_Handler::EmptyEventPhases() 
{
  if (p_phases) {
    while (!p_phases->empty()) {
      delete p_phases->back();
      p_phases->pop_back();
    }
  }
}  

void Event_Handler::PrintGenericEventStructure() 
{
  if (!msg_LevelIsInfo()) return;
  msg_Out()<<"----------------------------------------------------------"<<std::endl
	    <<"-- SHERPA generates events with the following structure --"<<std::endl
	    <<"----------------------------------------------------------"<<std::endl;
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg_Out()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  msg_Out()<<"---------------------------------------------------------"<<std::endl;
}

void Event_Handler::Reset(const bool sameevent)
{
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) 
    if (!sameevent || (*pit)->Type()!=eph::Perturbative ||
	(*pit)->Name().find("Signal_Processes")==std::string::npos) 
      (*pit)->CleanUp();
  if (!sameevent) {
    m_blobs.Clear();
    if (Particle::Counter()>m_lastparticlecounter || 
	Blob::Counter()>m_lastblobcounter) {
      msg_Error()<<METHOD<<"(): "<<Particle::Counter()
		 <<" particles and "<<Blob::Counter()
		 <<" blobs undeleted. Continuing."<<std::endl;
      m_lastparticlecounter=Particle::Counter();
      m_lastblobcounter=Blob::Counter();
    }

    Blob::Reset();
    Particle::Reset();
    Flow::ResetCounter();

    Blob *signal(new Blob());
    signal->SetType(btp::Signal_Process);
    signal->SetId();
    signal->SetStatus(blob_status::needs_signal);
    m_blobs.push_back(signal);
  }
  else {
    if (p_mehandler && p_mehandler->Weight()!=1.0) p_mehandler->SaveNumberOfTrials();
    Blob *signal(m_blobs.FindFirst(btp::Signal_Process));
    m_blobs.Clear(signal);
    signal->SetStatus(blob_status::internal_flag |
		      blob_status::needs_signal);
    Blob::Reset();
    Particle::Reset();
    Flow::ResetCounter();
  }
} 

bool Event_Handler::GenerateEvent(int mode, bool reset) 
{
  ATOOLS::ran.SaveStatus();
#ifdef USING__PYTHIA
  Lund_Interface::SaveStatus();
#endif
  PROFILE_LOCAL("Event_Handler::GenerateEvent");
  if (!rpa.gen.CheckTime()) {
    msg_Error()<<ATOOLS::om::bold
                     <<"\n\nEvent_Handler::GenerateEvent("<<mode<<"): "
                     <<ATOOLS::om::reset<<ATOOLS::om::red
                     <<"Timeout. Interrupt event generation."
                     <<ATOOLS::om::reset<<std::endl;
    kill(getpid(),SIGINT);
  }
  double weight = 1.;
  switch (mode) {
  case 0:
    if(reset) Reset();
    do {
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();) {
	if ((*pit)->Type()==eph::Analysis) {
	  ++pit;
	  continue;
	}
	msg_Debugging()<<METHOD<<"(): run '"<<(*pit)->Name()<<"'"<<std::endl;
	switch ((*pit)->Treat(&m_blobs,weight)) {
	case Return_Value::Nothing :
	  ++pit;
	  break;
	case Return_Value::Success : 
	  rvalue.IncCall((*pit)->Name());
	  pit=p_phases->begin();
	  break;
	case Return_Value::Error :
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncError((*pit)->Name());
	  return false;
	case Return_Value::Retry_Phase : 
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncRetryPhase((*pit)->Name());
	  pit=p_phases->begin();
	  break;
	case Return_Value::Retry_Event : 
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncRetryEvent((*pit)->Name());
	  Reset(true);
	  pit=p_phases->begin();
	  break;
	case Return_Value::New_Event : 
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncNewEvent((*pit)->Name());
	  Reset();
	  weight=1.;
	  pit=p_phases->begin();
	  break;
	default:
	  THROW(fatal_error,"Invalid return value");
	}
	if (weight==0.0 && rpa.gen.NumberOfDicedEvents()==
	    rpa.gen.NumberOfEvents()) return true;
      }
    } while (m_blobs.empty() || 
	     m_blobs.FindFirst(btp::Signal_Process)->NInP()==0);
    if (!m_blobs.FourMomentumConservation()) return false;
    if(p_mehandler) p_mehandler->ResetNumberOfTrials();
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) {
	switch ((*pit)->Treat(&m_blobs,weight)) {
	case Return_Value::Nothing :
	  break;
	case Return_Value::Success : 
	  rvalue.IncCall((*pit)->Name());
	  msg_Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "
			<<ATOOLS::om::reset
			<<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
			<<" yields "<<ATOOLS::om::bold<<true<<ATOOLS::om::reset<<std::endl;
	  break;
	case Return_Value::Error :
	  rvalue.IncCall((*pit)->Name());
	  rvalue.IncError((*pit)->Name());
	  return false;
	default:
	  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		     <<"  Unknown return value for 'Treat',"<<std::endl
		     <<"  Will continue and hope for the best."<<std::endl;
	  return false;
	}
      }
    }
    return true;
  case 9000:
  case 9001:
  case 9002: {
    if(reset) Reset();
    m_blobs.Clear();
    bool flag = false;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::External_MC) {
	bool result=(*pit)->Treat(&m_blobs,weight);
	if (result) flag = true;
      }
    }
    if (!flag) {
      msg_Error()<<"Error in Event_Handler::GenerateEvent"<<std::endl
		 <<"   Treat by External_MC failed, mode = "<<mode<<"."<<std::endl
		 <<"   Return 0 and hope for the best."<<std::endl;
    }
    if (m_blobs.empty()) {
      msg_Out()<<"Potential error in Event_Handler::GenerateEvent"<<std::endl
	       <<"   Empty bloblist would go into analysis !"<<std::endl;
      return false;
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }

    return true;
  }
  case 9999: {
    if(reset) Reset();
    bool flag(true);
    while (flag) {
      flag = false;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::Read_In) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  if (result) flag = true;
	}
      }
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }
    return true;
  }
  }
  return false;
} 

void Event_Handler::Finish() {
  msg_Info()<<"In Event_Handler::Finish : Summarizing the run may take some time."<<std::endl;
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
    (*pit)->Finish(std::string("Results"));
    (*pit)->CleanUp();
  }
  m_blobs.Clear();
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   After event : "<<Particle::Counter()<<" / "<<Blob::Counter()
	       <<" particles / blobs undeleted !"<<std::endl
	       <<"   Continue and hope for the best."<<std::endl;
    m_lastparticlecounter=Particle::Counter();
    m_lastblobcounter=Blob::Counter();
  }
  Blob::Reset();
  Particle::Reset();
  Flow::ResetCounter();
}
