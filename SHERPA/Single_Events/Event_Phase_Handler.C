#include "SHERPA/Single_Events/Event_Phase_Handler.H"

using namespace SHERPA;

std::ostream& SHERPA::operator<<(std::ostream& ostr, const eph::code ephc) {
  switch (ephc) {
  case eph::Unspecified:        return ostr<<"Unspecified       ";
  case eph::Perturbative:       return ostr<<"Perturbative      ";
  case eph::Hadronization:      return ostr<<"Hadronization     ";
  case eph::Analysis:           return ostr<<"Analysis          ";
  case eph::Read_In:            return ostr<<"Read_In           ";
  default:                      return ostr<<"Unknown           ";
  }
}

Event_Phase_Handler::Event_Phase_Handler() :
  m_type(eph::Unspecified), m_name(std::string("No Name")) { }

Event_Phase_Handler::Event_Phase_Handler(std::string _name) :
  m_type(eph::Unspecified), m_name(_name) { }


Event_Phase_Handler::~Event_Phase_Handler() { }

