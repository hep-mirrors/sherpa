//bof
//Version: 2 ADICIC++-0.0/2004/08/20

//Implementation of Dipole_Parameter.H.



#include <cassert>
#include <cstdlib>
#include <iostream>
#include "Recoil_Strategy.hpp"
#include "Evolution_Strategy.hpp"
#include "Dipole_Parameter.H"





using namespace std;
using namespace ADICIC;





//#include "..."





//=============================================================================



bool Dipole_Parameter::s_isalphasrun=false;
double Dipole_Parameter::s_alphasfix=0.0;
double Dipole_Parameter::s_k2tmin=0.0;
double Dipole_Parameter::s_k2tmax=0.0;

int Dipole_Parameter::s_restratqqbar=0;
int Dipole_Parameter::s_restratqg   =0;
int Dipole_Parameter::s_restratgqbar=0;
int Dipole_Parameter::s_restratgg   =0;

int Dipole_Parameter::s_chevolstrat=0;

const bool Dipole_Parameter::sf_start=Dipole_Parameter::ForceFirstInit();



//=============================================================================



Dipole_Parameter::Dipole_Parameter() {}

Dipole_Parameter::~Dipole_Parameter() {}



//=============================================================================



void Dipole_Parameter::Show() {    //Static.
  cout<<endl;
  cout<<"============================================"<<endl;
  cout<<"    Current ADICIC parameter adjustments"<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<"Enable possibility of running AlphaS = "<<s_isalphasrun<<"."<<endl;
  cout<<"The Value for fixed AlphaS  = "<<s_alphasfix<<"."<<endl;
  cout<<"Dipole shower cut-off scale = "<<s_k2tmin<<" GeV^2."<<endl;
  cout<<"Dipole shower maximum scale = "<<s_k2tmax<<" GeV^2."<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<"Recoil strategy for q-qbar dipoles = "<<s_restratqqbar<<"."<<endl;
  cout<<"Recoil strategy for q-g dipoles    = "<<s_restratqg<<"."<<endl;
  cout<<"Recoil strategy for g-qbar dipoles = "<<s_restratgqbar<<"."<<endl;
  cout<<"Recoil strategy for g-g dipoles    = "<<s_restratgg<<"."<<endl;
  cout<<"--------------------------------------------"<<endl;
  cout<<"Chain evolution strategy           = "<<s_chevolstrat<<"."<<endl;
  cout<<"============================================"<<endl;
}





const bool Dipole_Parameter::Check(const double cmsen) {    //Static.
  if(0.001<s_alphasfix && s_alphasfix<1.0); else {
    cerr<<"\nParameter out of range: value of fixed AlphaS!\n";
    assert(0.001<s_alphasfix && s_alphasfix<1.0);
  }
  if(0.0<s_k2tmin && s_k2tmin<10.0); else {
    cerr<<"\nParameter out of range: dipole shower cut-off scale!\n";
    assert(0.0<s_k2tmin && s_k2tmin<10.0);
  }
  if(s_k2tmin<s_k2tmax); else {
    cerr<<"\nParameter out of range: dipole shower maximum scale!\n";
    assert(s_k2tmin<s_k2tmax);
  }
  if(s_k2tmax<=cmsen*cmsen); else {
    cerr<<"\nParameter out of range: dipole shower maximum scale!\n";
    assert(s_k2tmax<=cmsen*cmsen);
  }
  return true;
}





const bool Dipole_Parameter::ForceFirstInit() {    //Static.
  static bool firsttime=true;
  if(firsttime==false) return false;
#ifdef DIPOLE_PARAMETER_OUTPUT
  cout<<"ADICIC::Dipole_Parameter::ForceFirstInit() is running.\n";
#endif
  firsttime=false;
#include "Dipole_Parameter.prm.cc"
  return Check(sqrt(s_k2tmax));
}



//=============================================================================



const bool Dipole_Parameter::Reset() {    //Static.
#include "Dipole_Parameter.prm.cc"
  return Check(sqrt(s_k2tmax));
}



//=============================================================================





//eof
