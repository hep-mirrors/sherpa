#include "Event_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"
#include "My_Limits.H"
#include <unistd.h>

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
  m_lastblobcounter(0)
{
  p_phases  = new Phase_List;
}

Event_Handler::~Event_Handler() 
{
  CleanUpEvent();
  EmptyEventPhases();
  
  if (p_phases)   { delete p_phases;   p_phases   = NULL; }
}

void Event_Handler::AddEventPhase(Event_Phase_Handler * phase) 
{
  eph::code type   = phase->Type();
  std::string name = phase->Name();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) { 
    if ((type==(*pit)->Type()) && (name==(*pit)->Name())) {
      msg.Out()<<"WARNING in Event_Handler::AddEventPhase("<<type<<":"<<name<<") "
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
  msg.Error()<<"Error in Event_Handler::GetEventPhase("<<i<<")"<<std::endl
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
  if (!msg.LevelIsInfo()) return;
  msg.Out()<<"----------------------------------------------------------"<<std::endl
	    <<"-- SHERPA generates events with the following structure --"<<std::endl
	    <<"----------------------------------------------------------"<<std::endl;
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg.Out()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  msg.Out()<<"---------------------------------------------------------"<<std::endl;
}

bool Event_Handler::GenerateEvent(int mode) 
{
  PROFILE_LOCAL("Event_Handler::GenerateEvent");
  if (!rpa.gen.CheckTime()) {
    ATOOLS::msg.Error()<<ATOOLS::om::bold
                     <<"\n\nEvent_Handler::GenerateEvent("<<mode<<"): "
                     <<ATOOLS::om::reset<<ATOOLS::om::red
                     <<"Timeout. Interrupt event generation."
                     <<ATOOLS::om::reset<<std::endl;
    kill(getpid(),SIGINT);
  }
  bool   flag   = true;
  double weight = 1.;
  Blob * hardblob;
  switch (mode) {
  case 0:
    do {
      CleanUpEvent();
      ATOOLS::ran.SaveStatus();
      hardblob = new Blob();
      hardblob->SetType(btp::Signal_Process);
      hardblob->SetStatus(-1);
      hardblob->SetId();
      hardblob->SetStatus(2);
      m_blobs.push_back(hardblob);
      while (flag) {
	flag = false;
	for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	  if ((*pit)->Type()==eph::Perturbative) {
	    bool result=(*pit)->Treat(&m_blobs,weight);
	    if (weight==0.0 &&
		rpa.gen.NumberOfDicedEvents()==rpa.gen.NumberOfEvents()) return true;
	    msg_Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "
			  <<ATOOLS::om::reset
			  <<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
			  <<" yields "<<ATOOLS::om::bold<<result<<ATOOLS::om::reset<<std::endl;
	    if (result) flag = true;
	  }
	}
      }
      flag=1;
      while (flag) {
	flag = false;
	for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	  if ((*pit)->Type()==eph::Hadronization) {
	    bool result=(*pit)->Treat(&m_blobs,weight);
	    msg_Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "<<ATOOLS::om::reset
				  <<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
				  <<" yields "<<ATOOLS::om::bold<<result<<ATOOLS::om::reset<<std::endl;
	    if (result) flag = true;
	  }
	}
      }
      flag=true;
    } while (m_blobs.empty());
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }
    return true;
  case 9000:
  case 9001:
  case 9002:
    CleanUpEvent();
    flag = false;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::External_MC) {
	bool result=(*pit)->Treat(&m_blobs,weight);
	if (result) flag = true;
      }
    }
    if (!flag) {
      msg.Error()<<"Error in Event_Handler::GenerateEvent"<<std::endl
		 <<"   Treat by External_MC failed, mode = "<<mode<<"."<<std::endl
		 <<"   Return 0 and hope for the best."<<std::endl;
    }
    if (m_blobs.empty()) {
      msg.Out()<<"Potential error in Event_Handler::GenerateEvent"<<std::endl
	       <<"   Empty bloblist would go into analysis !"<<std::endl;
      return false;
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }

    return true;
  case 9999:
    CleanUpEvent();
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
  return false;
} 


void Event_Handler::CleanUpEvent() 
{
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      (*pit)->CleanUp();
    }
  }
  m_blobs.Clear();
  if (Particle::Counter()>m_lastparticlecounter || 
      Blob::Counter()>m_lastblobcounter) {
    msg.Error()<<"Error in Event_Handler::CleanUpEvent()"<<std::endl
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




void Event_Handler::Finish() {
  msg_Info()<<"In Event_Handler::Finish : Summarizing the run may take some time."<<std::endl;
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) 
    (*pit)->Finish(std::string("Results"));
}

