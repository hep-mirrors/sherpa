#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace CFPSHOWER;
using namespace std;

ostream & CFPSHOWER::operator<<(ostream &s,const kernel_type::code & type) {
  switch (int(type)) {
  case 1: s<<"FF";break;
  case 2: s<<"FI";break;
  case 3: s<<"IF";break;
  case 4: s<<"II";break;
  case 0:
  default: s<<"unknown";break;
  }
  return s;
}

ostream & CFPSHOWER::operator<<(ostream &s,const Kernel_Info & info) {
  if (info.m_flavs.size()==0) s<<"  Kernel for undefined flavours: ";
  else {
    s<<"  Kernel for ["<<info.m_flavs[0].Bar()<<" -->";
    for (size_t i=1;i<info.m_flavs.size();i++) s<<" "<<info.m_flavs[i];
    s<<"]: ";
  }
  s<<info.m_type<<".\n";
  return s;
}

Kernel_Info::Kernel_Info(MODEL::Single_Vertex * vertex,
			 ATOOLS::Flavour_Vector & flavs,
			 kernel_type::code type,const bool & swapped) :
  p_alphaS(NULL), p_alpha(NULL),
  p_vertex(vertex), m_flavs(flavs), m_type(type), m_swapped(swapped) {
  m_flavs[0]=m_flavs[0].Bar();
  if (m_type==kernel_type::FF &&
      m_flavs[0].IsFermion() && 
      m_flavs[1].IsVector() &&
      m_flavs[2].IsFermion()) {
    swap(m_flavs[1], m_flavs[2]);
  }
}


const string Kernel_Info::SFName() const {
  string name = "";
  switch (m_type) {
  case kernel_type::FF:   name += "FF_"; break;
  case kernel_type::FI:   name += "FI_"; break;
  case kernel_type::IF:   name += "IF_"; break;
  case kernel_type::II:   name += "II_"; break;
  case kernel_type::none: 
  default:
    break;
  }
  for (size_t i=0;i<m_flavs.size();i++) {
    switch (m_flavs[i].IntSpin()) {
    case 2: name += "V"; break;
    case 1: name += "F"; break;
    case 0: name += "S"; break;
    }
  }
  return name;
}
  
const string Kernel_Info::GPName() const {
  string name = "";
  for (size_t i=0;i<m_flavs.size();i++) {
    switch (abs(m_flavs[i].StrongCharge())) {
    case 8: name += "G"; break;
    case 3: name += "Q"; break;
    }
  }
  return name;
}
