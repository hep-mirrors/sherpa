#include "Event_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"

#include "Random.H"

#ifdef PROFILE__Event_Handler
#include "prof.hh"
#else 
#define PROFILE_HERE {}
#define PROFILE_LOCAL(LOCALNAME) {}
#endif

using namespace SHERPA;
using namespace ATOOLS;


Event_Handler::Event_Handler() 
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
      msg.Events()<<"Event_Handler::AddEventPhase("<<type<<":"<<name<<") "
		  <<"already included."<<std::endl;
      return;
    }
  }
  msg.Events()<<"Event_Handler::AddEventPhase("<<type<<":"<<name<<")."<<std::endl;
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
  msg.Events()<<"----------------------------------------------------------"<<std::endl
	      <<"-- SHERPA generates events with the following structure --"<<std::endl
	      <<"----------------------------------------------------------"<<std::endl;
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg.Events()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
  msg.Events()<<"---------------------------------------------------------"<<std::endl;
}

bool Event_Handler::GenerateEvent(int mode) 
{
  PROFILE_LOCAL("Event_Handler::GenerateEvent");

  CleanUpEvent();

  bool flag     = 1;
  double weight = 1.;
  Blob * hardblob;
  switch (mode) {
  case 0:
    hardblob = new Blob();
    hardblob->SetType(btp::Signal_Process);
    hardblob->SetStatus(-1);
    hardblob->SetId(0);
    hardblob->SetStatus(2);
    m_blobs.push_back(hardblob);
    while (flag) {
      flag = 0;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::Perturbative) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  ATOOLS::msg.Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "<<ATOOLS::om::reset
				<<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
				<<" yields "<<ATOOLS::om::bold<<result<<ATOOLS::om::reset<<std::endl;
	  if (result) flag = 1;
	}
      }
    }
    flag=1;
    while (flag) {
      flag = 0;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::Hadronization) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  ATOOLS::msg.Tracking()<<ATOOLS::om::blue<<"Event_Handler::GenerateEvent("<<mode<<"): "<<ATOOLS::om::reset
				<<"Event phase "<<ATOOLS::om::bold<<(*pit)->Name()<<ATOOLS::om::reset
				<<" yields "<<ATOOLS::om::bold<<result<<ATOOLS::om::reset<<std::endl;
	  if (result) flag = 1;
	}
      }
    }
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }
    PrintBlobs();
    return 1;
  case 9000:
    while (flag) {
      flag = 0;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::External_MC) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  if (result) flag = 1;
	}
      }
    }
    
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }

    return 1;
  case 9999:
    while (flag) {
      flag = 0;
      for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
	if ((*pit)->Type()==eph::Read_In) {
	  bool result=(*pit)->Treat(&m_blobs,weight);
	  if (result) flag = 1;
	}
      }
    }      
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==eph::Analysis) (*pit)->Treat(&m_blobs,weight);
    }
    return 1;
  }
  return 0;
} 


void Event_Handler::CleanUpEvent() 
{
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      (*pit)->CleanUp();
    }
  }
  if (!m_blobs.empty()) {
    for (Blob_Iterator blit=m_blobs.begin();blit!=m_blobs.end();++blit) delete (*blit);
    m_blobs.clear();
  }
  Flow::ResetCounter();
  Particle::Reset();
}


void Event_Handler::PrintEvent(int mode) {
  switch (mode) {
    case 1:  PrintBlobs();   return;
  }
}

void Event_Handler::PrintBlobs() {
  if (!msg.LevelIsEvents()) return;
  msg.Out()<<"  -------------------------------------------------  "<<std::endl;
  if (!m_blobs.empty()) {
    for (Blob_Iterator blit=m_blobs.begin();blit!=m_blobs.end();++blit) {
      msg.Out()<<(*blit)<<std::endl;
    }
  }
  else {
    msg.Out()<<"  ***** Empty Event *****  "<<std::endl;
  }
  msg.Out()<<"  -------------------------------------------------  "<<std::endl;
}




void Event_Handler::Finish() {
  std::cout<<"In SummarizeRun()"<<std::endl;
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) 
    (*pit)->Finish(std::string("Results"));
}

