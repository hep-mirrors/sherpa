//bof
//Version: 1 ADICIC++-0.0/2004/03/03

//Implementation of Dipole_Flavour.H.



#include "Dipole_Flavour.H"





using namespace ATOOLS;
using namespace ADICIC;





const bool Dipole_Flavour_Init::status=Dipole_Flavour_Init::Now();



//gluon-type-like tags: even
//quark-like tags     : odd
//antiquark-like tags : -1 * quark-like tag



const short Dipole_Gluon_Base::GNtag=0;
const Flavour Dipole_Gluon_Base::GN(kf::gluon);



const short Dipole_Quark_Base::QUtag=1;
const Flavour Dipole_Quark_Base::DQ(kf::d);
const Flavour Dipole_Quark_Base::UQ(kf::u);
const Flavour Dipole_Quark_Base::SQ(kf::s);
const Flavour Dipole_Quark_Base::CQ(kf::c);
const Flavour Dipole_Quark_Base::BQ(kf::b);
const Flavour Dipole_Quark_Base::TQ(kf::t);



const short Dipole_Antiquark_Base::AQtag=-1;
const Flavour Dipole_Antiquark_Base::DA(kf::d,1);
const Flavour Dipole_Antiquark_Base::UA(kf::u,1);
const Flavour Dipole_Antiquark_Base::SA(kf::s,1);
const Flavour Dipole_Antiquark_Base::CA(kf::c,1);
const Flavour Dipole_Antiquark_Base::BA(kf::b,1);
const Flavour Dipole_Antiquark_Base::TA(kf::t,1);





//eof











