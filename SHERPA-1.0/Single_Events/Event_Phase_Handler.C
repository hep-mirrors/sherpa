#include "Event_Phase_Handler.H"
#include "Message.H"

using namespace SHERPA;
using namespace ATOOLS;

Event_Phase_Handler::Event_Phase_Handler() :
  m_type(std::string("Unspecified")), m_name(std::string("No Name")) { }

Event_Phase_Handler::Event_Phase_Handler(std::string _name) :
  m_type(std::string("Unspecified")), m_name(_name) { }


Event_Phase_Handler::~Event_Phase_Handler() { }

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
