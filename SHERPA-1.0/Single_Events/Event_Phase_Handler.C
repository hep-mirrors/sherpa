#include "Event_Phase_Handler.H"

using namespace SHERPA;

Event_Phase_Handler::Event_Phase_Handler() :
  m_type(std::string("Unspecified")), m_name(std::string("No Name")) { }

Event_Phase_Handler::Event_Phase_Handler(std::string _name) :
  m_type(std::string("Unspecified")), m_name(_name) { }


Event_Phase_Handler::~Event_Phase_Handler() { }

