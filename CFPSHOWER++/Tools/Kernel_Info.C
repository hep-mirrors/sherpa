#include "CFPSHOWER++/Tools/Kernel_Info.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>
#include <iostream>

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

ostream & CFPSHOWER::operator<<(ostream &s,const log_type::code & type) {
  switch (int(type)) {
  case 1: s<<"soft";break;
  case 2: s<<"coll";break;
  case 0: s<<"none";break;
  default: s<<"unknown";break;
  }
  return s;
}

ostream & CFPSHOWER::operator<<(ostream &s,const kin_type::code & type) {
  switch (int(type)) {
  case 2: s<<"Catani-Seymour";break;
  case 1: s<<"Alaric";break;
  default: s<<"unknown";break;
  }
  return s;
}

ostream & CFPSHOWER::operator<<(ostream &s,const muR_scheme::code & scheme) {
  switch (int(scheme)) {
  case 11: s<<"muR^2 = KT^2 (M^2 for g->qq)";break;
  case 10: s<<"muR^2 = KT^2";break;
  case 1:  s<<"muR^2 = t (M^2 for g->qq)";break;
  case 0:  s<<"muR^2 = t";break;
  default: s<<"unknown";break;
  }
  return s;
}


ostream & CFPSHOWER::operator<<(ostream &s,const Kernel_Info & info) {
  if (info.m_flavs.size()==0) s<<"  Kernel for undefined flavours\n";
  else {
    s<<"  Kernel ["<<info.m_split<<" -->";
    for (size_t i=0;i<info.m_flavs.size();i++)
      s<<" "<<info.m_flavs[i];
    s<<"] [";
    for (size_t i=0;i<info.m_flavs.size();i++)
      s<<" "<<info.m_tagsequence[i];
    s<<"]\n";
  }
  return s;
}

Kernel_Info::Kernel_Info(log_type::code logtype,
			 kin_type::code kintype,
			 kernel_type::code type) :
  m_logtype(logtype), m_kintype(kintype), m_type(type) {}

Kernel_Info::Kernel_Info(ATOOLS::Flavour & split,ATOOLS::Flavour_Vector & flavs,
			 const vector<size_t> & tagsequence,
			 log_type::code logtype,
			 kin_type::code kintype,
			 muR_scheme::code muRscheme,
			 kernel_type::code type) :
  m_split(split), m_flavs(flavs), m_tagsequence(tagsequence),
  m_logtype(logtype), m_kintype(kintype), m_muRscheme(muRscheme), m_type(type),
  p_alphaS(NULL), p_alpha(NULL),
  m_kfactor(0), m_softcorr(0), m_endpoint(0),
  m_asfactor(1.), m_muR2factor(1.)
{
  AdjustFlavours();
}

void Kernel_Info::AdjustFlavours() {
  //msg_Out()<<METHOD<<" for "<<m_split<<" -> ";
  //for (size_t i=0;i<m_flavs.size();i++) msg_Out()<<m_flavs[i];
  //msg_Out()<<"\n";
  if (m_type==kernel_type::FF || m_type==kernel_type::FI) {
    if (m_flavs[0]!=m_split) { 
      for (size_t i=1;i<m_flavs.size();i++) {
	if (m_flavs[i]==m_split) {
	  swap(m_flavs[0], m_flavs[i]);
	  swap(m_tagsequence[0], m_tagsequence[i]);
	}
      }
    }
    size_t off = (m_split.IsVector() && m_flavs.size()==2) ? 0 : 1;
    if (!(off>0 && m_flavs.size()<3)) {
      size_t pos = (m_split.IsFermion() && m_split.IsAnti()) ? 0 : 1;
      size_t i = off+(1-pos), j = off+pos;
      if (m_flavs[i].IsFermion() && m_flavs[j].IsFermion() && m_flavs[j].IsAnti()) {
	swap(m_flavs[i], m_flavs[j]);
	swap(m_tagsequence[i], m_tagsequence[j]);
      }
    }
  }
}

void Kernel_Info::AdjustFlavours12() {
  if (m_type==kernel_type::FF || m_type==kernel_type::FI) {
    // make sure we have always F -> FV in terms of flavour sequence
    if (m_split.IsFermion() && m_flavs[0].IsVector() && m_flavs[1].IsFermion()) {
      swap(m_flavs[0], m_flavs[1]);
    }
    // make sure we have always V -> Fbar F in terms of flavour sequence
    if (m_split.IsVector() && m_flavs[0].IsFermion() &&
	m_flavs[1].IsFermion() && m_flavs[1].IsAnti()) {
      swap(m_flavs[0], m_flavs[1]);
    }
  }
  else if (m_type==kernel_type::IF) {
    if (m_tagsequence[0]==2) {
      m_split = m_split.Bar();
      for (size_t i=0;i<m_flavs.size();i++) m_flavs[i] = m_flavs[i].Bar(); 
    }
  }
}

const string Kernel_Info::KinName() const {
  string name = "Kinematics_";
  switch (m_type) {
  case kernel_type::FF: name += "FF"; break;
  case kernel_type::FI: name += "FI"; break;
  case kernel_type::IF: name += "IF"; break;
  case kernel_type::II: name += "II"; break;
  case kernel_type::none: 
  default:
    break;
  }
  name += "_";
  switch (m_kintype) {
  case kin_type::CataniSeymour:
    name += "Catani-Seymour"; break;
  case kin_type::Alaric:
    name += "Alaric";
    break;
  case kernel_type::none:
    name += "None"; break;
  default:
    break;
  }
  return name;
}

const string Kernel_Info::Name() const {
  string name = "Splitting_";
  switch (m_type) {
  case kernel_type::FF: name += "FF"; break;
  case kernel_type::FI: name += "FI"; break;
  case kernel_type::IF: name += "IF"; break;
  case kernel_type::II: name += "II"; break;
  case kernel_type::none: 
  default:
    break;
  }
  name += to_string(m_flavs.size());
  return name;
}
  
const string Kernel_Info::SFName() const {
  string name = "";
  switch (m_type) {
  case kernel_type::FF: name += "FF_"; break;
  case kernel_type::FI: name += "FI_"; break;
  case kernel_type::IF: name += "IF_"; break;
  case kernel_type::II: name += "II_"; break;
  case kernel_type::none: 
  default:
    break;
  }
  switch (m_logtype) {
  case log_type::soft: name += "Soft"; break;
  case log_type::coll:
    name += "Coll_";
    switch (m_split.IntSpin()) {
    case 2: name += "V"; break;
    case 1: name += "F"; break;
    case 0: name += "S"; break;
    }
    for (size_t i=0;i<m_flavs.size();i++) {
      switch (m_flavs[i].IntSpin()) {
      case 2: name += "V"; break;
      case 1: name += "F"; break;
      case 0: name += "S"; break;
      }
    }
    break;
  case log_type::none: 
  default:
    break;
  }
  return name;
}
  
const string Kernel_Info::GPName() const {
  string name = "";
  switch (abs(m_split.StrongCharge())) {
  case 8: name += "G"; break;
  case 3: name += "Q"; break;
  }
  for (size_t i=0;i<m_flavs.size();i++) {
    switch (abs(m_flavs[i].StrongCharge())) {
    case 8: name += "G"; break;
    case 3: name += "Q"; break;
    }
  }
  return name;
}
