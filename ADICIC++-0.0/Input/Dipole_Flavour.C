//bof
//Version: 1 ADICIC++-0.0/2004/07/05

//Implementation of Dipole_Flavour.H.



#include "Dipole_Flavour.H"





using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



bool Dipole_Flavour_Init::s_status=false;


const bool Dipole_Flavour_Init::DoIt(bool level) {    //Static.

  if(s_status) return false;
  if(level) ATOOLS::ParticleInit("../TestIt/data");

  static const Flavour gn(kf::gluon); Dipole_Gluon_Base::GN=&gn;

  static const Flavour dq(kf::d); Dipole_Quark_Base::DQ=&dq;
  static const Flavour uq(kf::u); Dipole_Quark_Base::UQ=&uq;
  static const Flavour sq(kf::s); Dipole_Quark_Base::SQ=&sq;
  static const Flavour cq(kf::c); Dipole_Quark_Base::CQ=&cq;
  static const Flavour bq(kf::b); Dipole_Quark_Base::BQ=&bq;
  static const Flavour tq(kf::t); Dipole_Quark_Base::TQ=&tq;

  static const Flavour da(kf::d,1); Dipole_Antiquark_Base::DA=&da;
  static const Flavour ua(kf::u,1); Dipole_Antiquark_Base::UA=&ua;
  static const Flavour sa(kf::s,1); Dipole_Antiquark_Base::SA=&sa;
  static const Flavour ca(kf::c,1); Dipole_Antiquark_Base::CA=&ca;
  static const Flavour ba(kf::b,1); Dipole_Antiquark_Base::BA=&ba;
  static const Flavour ta(kf::t,1); Dipole_Antiquark_Base::TA=&ta;

  s_status=true;
  return true;

}



//=============================================================================



//gluon-type-like tags: even
//quark-like tags     : odd
//antiquark-like tags : -1 * quark-like tag



const short Dipole_Gluon_Base::GNtag=0;
const Flavour* Dipole_Gluon_Base::GN(NULL);



const short Dipole_Quark_Base::QUtag=1;
const Flavour* Dipole_Quark_Base::DQ(NULL);
const Flavour* Dipole_Quark_Base::UQ(NULL);
const Flavour* Dipole_Quark_Base::SQ(NULL);
const Flavour* Dipole_Quark_Base::CQ(NULL);
const Flavour* Dipole_Quark_Base::BQ(NULL);
const Flavour* Dipole_Quark_Base::TQ(NULL);



const short Dipole_Antiquark_Base::AQtag=-1;
const Flavour* Dipole_Antiquark_Base::DA(NULL);
const Flavour* Dipole_Antiquark_Base::UA(NULL);
const Flavour* Dipole_Antiquark_Base::SA(NULL);
const Flavour* Dipole_Antiquark_Base::CA(NULL);
const Flavour* Dipole_Antiquark_Base::BA(NULL);
const Flavour* Dipole_Antiquark_Base::TA(NULL);



//=============================================================================



Dipole_Gluon_G::Dipole_Gluon_G() {}
Dipole_Gluon_G::~Dipole_Gluon_G() {}


Dipole_Quark_D::Dipole_Quark_D() {}
Dipole_Quark_D::~Dipole_Quark_D() {}
Dipole_Quark_U::Dipole_Quark_U() {}
Dipole_Quark_U::~Dipole_Quark_U() {}
Dipole_Quark_S::Dipole_Quark_S() {}
Dipole_Quark_S::~Dipole_Quark_S() {}
Dipole_Quark_C::Dipole_Quark_C() {}
Dipole_Quark_C::~Dipole_Quark_C() {}
Dipole_Quark_B::Dipole_Quark_B() {}
Dipole_Quark_B::~Dipole_Quark_B() {}
Dipole_Quark_T::Dipole_Quark_T() {}
Dipole_Quark_T::~Dipole_Quark_T() {}


Dipole_Antiquark_D::Dipole_Antiquark_D() {}
Dipole_Antiquark_D::~Dipole_Antiquark_D() {}
Dipole_Antiquark_U::Dipole_Antiquark_U() {}
Dipole_Antiquark_U::~Dipole_Antiquark_U() {}
Dipole_Antiquark_S::Dipole_Antiquark_S() {}
Dipole_Antiquark_S::~Dipole_Antiquark_S() {}
Dipole_Antiquark_C::Dipole_Antiquark_C() {}
Dipole_Antiquark_C::~Dipole_Antiquark_C() {}
Dipole_Antiquark_B::Dipole_Antiquark_B() {}
Dipole_Antiquark_B::~Dipole_Antiquark_B() {}
Dipole_Antiquark_T::Dipole_Antiquark_T() {}
Dipole_Antiquark_T::~Dipole_Antiquark_T() {}



//=============================================================================





//eof
