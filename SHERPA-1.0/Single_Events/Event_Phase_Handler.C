#include "Event_Phase_Handler.H"
#include "Message.H"

using namespace SHERPA;
using namespace AORGTOOLS;

Event_Phase_Handler::Event_Phase_Handler()
{
  m_type = std::string("Unspecified");
  m_name = std::string("No Name");
}


Event_Phase_Handler::~Event_Phase_Handler() {
  EmptyMyLists();
}

void Event_Phase_Handler::EmptyMyLists() {
  msg.Debugging()<<"Try to delete blobs in "<<m_name<<" : "<<m_myblobs.size()<<endl;
  while (!m_myblobs.empty()) {
    delete m_myblobs.back();
    m_myblobs.pop_back();
  }
  /*
    msg.Debugging()<<"Try to delete partons in "<<m_name<<" : "<<m_mypartons.size()<<endl;
    while (!m_mypartons.empty()) {
    delete m_mypartons.back();
    m_mypartons.pop_back();
    }
  */
}


int Event_Phase_Handler::NumberOfBlobs() {
  return m_myblobs.size(); 
}

APHYTOOLS::Blob_List * Event_Phase_Handler::GetBlobs() { 
  return &m_myblobs; 
}

APHYTOOLS::Blob * Event_Phase_Handler::GetBlob(int i) {
  //if ((i>-1) && (i<m_myblobs.size())) return m_myblobs[i];
  msg.Error()<<"Error in Event_Phase_Handler::Blob("<<i<<") : "
	     <<"Out of bounds : "<<m_myblobs.size()<<endl;
  return NULL;
} 

void Event_Phase_Handler::AddBlob(APHYTOOLS::Blob * _blob) {
  m_myblobs.push_back(_blob);
}

int Event_Phase_Handler::NumberOfPartons() {
  //return m_mypartons.size(); 
}

APHYTOOLS::Parton_List * Event_Phase_Handler::GetPartons() { 
  //return &m_mypartons; 
}

APHYTOOLS::Parton * Event_Phase_Handler::GetParton(int i) {
  //if ((i>-1) && (i<m_mypartons.size())) return m_mypartons[i];
  msg.Error()<<"Error in Event_Phase_Handler::Parton("<<i<<") : "
	     <<"Out of bounds : "<<m_mypartons.size()<<endl;
  return NULL;
} 

void Event_Phase_Handler::AddParton(APHYTOOLS::Parton * _parton) {
  m_mypartons.push_back(_parton);
}




std::string Event_Phase_Handler::Name() {
  return m_name;
}

void Event_Phase_Handler::SetName(std::string _name) {
  m_name = _name;
}

std::string Event_Phase_Handler::Type() {
  return m_type;
}

void Event_Phase_Handler::SetType(std::string _type) {
  m_type = _type;
}
