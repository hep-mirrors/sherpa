//bof
//Version: 3 ADICIC++-0.0/2005/07/23

//Implementation of Sudakov_Group.H.



#include "Sudakov_Group.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





#include "Sudakov_Group.tpt.cc"





//=============================================================================



template class Sudakov_Group<Dipole::qqbar>;
template class Sudakov_Group<Dipole::qg>;
template class Sudakov_Group<Dipole::gqbar>;
template class Sudakov_Group<Dipole::gg>;



//=============================================================================



template<>
const short Sudakov<Dipole::qqbar,Radiation::gluon>::s_x1pow=2;
template<>
const short Sudakov<Dipole::qqbar,Radiation::gluon>::s_x3pow=2;
template<>
const double Sudakov<Dipole::qqbar,Radiation::gluon>::s_colfac=1.5*M_PI;

template<>
const short Sudakov<Dipole::qg,Radiation::gluon>::s_x1pow=2;
template<>
const short Sudakov<Dipole::qg,Radiation::gluon>::s_x3pow=3;
template<>
const double Sudakov<Dipole::qg,Radiation::gluon>::s_colfac=4.0*M_PI/3.0;

template<>
const short Sudakov<Dipole::gqbar,Radiation::gluon>::s_x1pow=3;
template<>
const short Sudakov<Dipole::gqbar,Radiation::gluon>::s_x3pow=2;
template<>
const double Sudakov<Dipole::gqbar,Radiation::gluon>::s_colfac=4.0*M_PI/3.0;

template<>
const short Sudakov<Dipole::gg,Radiation::gluon>::s_x1pow=3;
template<>
const short Sudakov<Dipole::gg,Radiation::gluon>::s_x3pow=3;
template<>
const double Sudakov<Dipole::gg,Radiation::gluon>::s_colfac=4.0*M_PI/3.0;



template<>
const double Sudakov<Dipole::qg,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double Sudakov<Dipole::gqbar,Radiation::quark>::s_colfac=8.0*M_PI;
template<>
const double Sudakov<Dipole::gg,Radiation::quark>::s_colfac=4.0*M_PI;



template class Sudakov<Dipole::qqbar,Radiation::gluon>;
template class Sudakov<Dipole::qg,Radiation::gluon>;
template class Sudakov<Dipole::gqbar,Radiation::gluon>;
template class Sudakov<Dipole::gg,Radiation::gluon>;

template class Sudakov<Dipole::qg,Radiation::quark>;
template class Sudakov<Dipole::gqbar,Radiation::quark>;
template class Sudakov<Dipole::gg,Radiation::quark>;



//=============================================================================





//eof
