#include "Event_Handler.H"
#include "Message.H"
#include "Run_Parameter.H"

#include "Random.H"

using namespace SHERPA;
using namespace ATOOLS;


Event_Handler::Event_Handler() :
  p_analysis(NULL)
{
  p_phases  = new Phase_List;
}

Event_Handler::~Event_Handler() 
{
  CleanUpEvent();
  EmptyEventPhases();
  
  if (p_phases)   { delete p_phases;   p_phases   = NULL; }
  if (p_analysis) { delete p_analysis; p_analysis = NULL; }
}

void Event_Handler::AddEventPhase(Event_Phase_Handler * _phase) 
{
  std::string _type = _phase->Type();
  std::string _name = _phase->Name();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) { 
    if ((_type==(*pit)->Type()) && (_name==(*pit)->Name())) {
      msg.Events()<<"Event_Handler::AddEventPhase("<<_type<<":"<<_name<<") "
		  <<"already included."<<std::endl;
      return;
    }
  }
  msg.Events()<<"Event_Handler::AddEventPhase("<<_type<<":"<<_name<<")."<<std::endl;
  p_phases->push_back(_phase);
}

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
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg.Events()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<std::endl;
    }
  }
}

bool Event_Handler::GenerateEvent() 
{
  CleanUpEvent();
  Blob * hardblob = new Blob();
  hardblob->SetType(std::string("Signal Process : "));
  hardblob->SetId(0);
  m_blobs.push_back(hardblob);

  bool flag = 1;
  double weight=1.;
  while (flag) {
    flag = 0;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==std::string("Perturbative")) {
	bool result=(*pit)->Treat(&m_blobs,weight);
 	ATOOLS::msg.Tracking()<<(*pit)->Name()<<" yields "<<result<<std::endl;
 	if (result) flag = 1;
      }
    }
  }

  // Maybe a boost phase here ???

  if (flag==0) flag=1;
  while (flag) {
    flag = 0;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==std::string("Hadronization")) {
	if ((*pit)->Treat(&m_blobs,weight)) flag = 1;
      }
    }
  }

  if (rpa.gen.Events()) PrintBlobs();
  return 1;
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
}


void Event_Handler::PrintEvent(int mode) {
  switch (mode) {
    case 1:  PrintBlobs();   return;
  }
}

void Event_Handler::PrintBlobs() {
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


void Event_Handler::PerformAnalysis() {
  //if (p_analysis) p_analysis->DoAnalysis(&m_partons);
}

void Event_Handler::SetAnalysis(Sample_Analysis * _analysis) { _analysis = p_analysis; }
