#include "Event_Handler.H"
#include "Message.H"

using namespace SHERPA;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;


Event_Handler::Event_Handler() 
{
  p_phases  = new Phase_List;
}

Event_Handler::~Event_Handler() 
{
  cout<<"in  Event_Handler::~Event_Handler() "<<endl;

  CleanUpEvent();
  EmptyEventPhases();
  
  if (p_phases)  { delete p_phases;  p_phases  = NULL; }
  cout<<"out Event_Handler::~Event_Handler() "<<endl;
}

void Event_Handler::AddEventPhase(Event_Phase_Handler * _phase) 
{
  std::string _type = _phase->Type();
  std::string _name = _phase->Name();
  for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) { 
    if ((_type==(*pit)->Type()) && (_name==(*pit)->Name())) {
      msg.Events()<<"Event_Handler::AddEventPhase("<<_type<<":"<<_name<<") "
		  <<"already included."<<endl;
      return;
    }
  }
  msg.Events()<<"Event_Handler::AddEventPhase("<<_type<<":"<<_name<<")."<<endl;
  p_phases->push_back(_phase);
}

void Event_Handler::EmptyEventPhases() 
{
  if (p_phases) {
    while (!p_phases->empty()) {
      cout<<" deleting "<<p_phases->back()->Name()<<endl;
      delete p_phases->back();
      p_phases->pop_back();
    }
  }
}  

void Event_Handler::PrintGenericEventStructure() 
{
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      msg.Debugging()<<(*pit)->Type()<<" : "<<(*pit)->Name()<<endl;
    }
  }
}

bool Event_Handler::GenerateEvent() 
{
  CleanUpEvent();
  Blob * hardblob = new Blob();
  hardblob->SetType(string("Signal Process : "));
  m_blobs.push_back(hardblob);

  bool flag = 1;
  msg.Debugging()<<"#################################################################"<<endl
		 <<"#################################################################"<<endl;
  while (flag) {
    flag = 0;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==string("Perturbative")) {
	msg.Debugging()<<"#############"<<endl
		       <<"   Try "<<(*pit)->Name()<<" : "<<m_blobs.size()<<endl;
	if ((*pit)->Treat(&m_blobs)) flag = 1;
      }
    }
  }

  if (flag==0) flag=1;
  while (flag) {
    flag = 0;
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      if ((*pit)->Type()==string("Hadronization")) {
	msg.Debugging()<<"#############"
		       <<"   Try "<<(*pit)->Name()<<" : "<<m_blobs.size()<<endl;
	if ((*pit)->Treat(&m_blobs)) flag = 1;
      }
    }
  }

  PrintBlobs();
  return 1;
}

void Event_Handler::CleanUpEvent() 
{
  if (!p_phases->empty()) {
    for (Phase_Iterator pit=p_phases->begin();pit!=p_phases->end();++pit) {
      (*pit)->CleanUp();
      msg.Debugging()<<"Cleaned up in "<<(*pit)->Type()<<" : "<<(*pit)->Name()<<endl;
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
    default: PrintPartons(); return;
  }
}

void Event_Handler::PrintBlobs() {
  msg.Out()<<"  -------------------------------------------------  "<<endl;
  if (!m_blobs.empty()) {
    for (Blob_Iterator blit=m_blobs.begin();blit!=m_blobs.end();++blit) {
      msg.Out()<<(*blit)<<endl;
    }
  }
  else {
    msg.Out()<<"  ***** Empty Event *****  "<<endl;
  }
  msg.Out()<<"  -------------------------------------------------  "<<endl;
}

void Event_Handler::PrintPartons() {
  /*
    msg.Out()<<"  -------------------------------------------------  "<<endl;
    if (!m_partons.empty()) {
    for (Parton_Iterator plit=m_partons.begin();plit!=m_partons.end();++plit) {
    msg.Out()<<(*plit)<<endl;
    }
    }
    else {
    msg.Out()<<"  ***** Empty Event *****  "<<endl;
    }
    msg.Out()<<"  -------------------------------------------------  "<<endl;
  */
}


void Event_Handler::PerformAnalysis() {
  p_analysis->DoAnalysis(&m_partons);
}

void Event_Handler::SetAnalysis(Sample_Analysis * _analysis) { _analysis = p_analysis; }
