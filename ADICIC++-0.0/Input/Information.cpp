//bof
//Version: 3 ADICIC++-0.0/2005/08/18

//Implementation of Information.hpp.



#include "Information.hpp"





//using namespace





ADICIC::Dipole_Flavour_Info::Dipole_Flavour_Info() {
  for(int i=0; i<22; ++i) gluon.pkf[i]=NULL;
  gluon.pkf[21]=&gluon.g;
  //for(int i=0; i<7; ++i) quark.pkf[i]=NULL;
  quark.pkf[0]=NULL;
  quark.pkf[1]=&quark.d;
  quark.pkf[2]=&quark.u;
  quark.pkf[3]=&quark.s;
  quark.pkf[4]=&quark.c;
  quark.pkf[5]=&quark.b;
  quark.pkf[6]=&quark.t;
  //for(int i=0; i<7; ++i) antiq.pkf[i]=NULL;
  antiq.pkf[0]=NULL;
  antiq.pkf[1]=&antiq.d;
  antiq.pkf[2]=&antiq.u;
  antiq.pkf[3]=&antiq.s;
  antiq.pkf[4]=&antiq.c;
  antiq.pkf[5]=&antiq.b;
  antiq.pkf[6]=&antiq.t;
}


ADICIC::Dipole_Flavour_Info::~Dipole_Flavour_Info() {}


const ADICIC::Dipole_Flavour_Info ADICIC::info;





//eof
