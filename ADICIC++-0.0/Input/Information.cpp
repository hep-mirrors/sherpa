//bof
//Version: ADICIC++-0.0/2004/09/01

//Implementation of Information.hpp.



#include "Information.hpp"





//using namespace





ADICIC::Dipole_Flavour_Info::Dipole_Flavour_Info() {}
ADICIC::Dipole_Flavour_Info::~Dipole_Flavour_Info() {}


const ADICIC::Dipole_Flavour_Info ADICIC::info;





ADICIC::Dipquarkbox::Dipquarkbox()
  : m_kfqkbase() {
  using namespace ATOOLS;
  m_kfqkbase[kf::d]=&info.quark.d;
  m_kfqkbase[kf::u]=&info.quark.u;
  m_kfqkbase[kf::s]=&info.quark.s;
  m_kfqkbase[kf::c]=&info.quark.c;
  m_kfqkbase[kf::b]=&info.quark.b;
  m_kfqkbase[kf::t]=&info.quark.t;
}

ADICIC::Dipantiqbox::Dipantiqbox()
  : m_kfaqbase() {
  using namespace ATOOLS;
  m_kfaqbase[kf::d]=&info.antiq.d;
  m_kfaqbase[kf::u]=&info.antiq.u;
  m_kfaqbase[kf::s]=&info.antiq.s;
  m_kfaqbase[kf::c]=&info.antiq.c;
  m_kfaqbase[kf::b]=&info.antiq.b;
  m_kfaqbase[kf::t]=&info.antiq.t;
}





ADICIC::Dipole_Flavour_Interface::Dipole_Flavour_Interface()
  : m_qk(), m_aq(), quark(m_qk), antiq(m_aq) {}
ADICIC::Dipole_Flavour_Interface::~Dipole_Flavour_Interface() {}


const ADICIC::Dipole_Flavour_Interface ADICIC::interface;





//eof
