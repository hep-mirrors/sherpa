//bof
//Version: 2 ADICIC++-0.0/2004/08/11

//Implementation of Paraminit.H.



#include "Run_Parameter.H"
#include "Dipole_Parameter.H"
#include "Sudakov_Calculator.H"
#include "Dipole_Handler.H"
#include "Chain_Handler.H"
#include "Paraminit.H"





using namespace std;
using namespace ATOOLS;
using namespace ADICIC;





//=============================================================================



unsigned Dipole_Parameter_Init::s_status=0;



//=============================================================================



const bool Dipole_Parameter_Init::DoIt() {    //Static.

  if(s_status) return false;

  ++s_status;

#ifdef PARAMINIT_OUTPUT
  Dipole_Parameter::Show();
#endif

  //###########################################################################
  //AlphaS treatment flag
  Dipole_Parameter::s_isalphasrun=true;//false;//true;
  //coupling
  Dipole_Parameter::s_alphasfix=0.15;//0.12;
  //GeV^2
  Dipole_Parameter::s_k2tmin=1.0;
  Dipole_Parameter::s_k2tmax=8100.0;
  //recoil strategies (compare with Recoil_Strategy.hpp)
  Dipole_Parameter::s_restratqqbar=2;    //7;//2;
  Dipole_Parameter::s_restratqg   =3;    //6;//3;
  Dipole_Parameter::s_restratgqbar=1;    //2;//1;
  Dipole_Parameter::s_restratgg   =4;//7;//5;//4;
  //chain evolution strategy (compare with Evolution_Strategy.hpp)
  Dipole_Parameter::s_chevolstrat =1;//2;//1;
  //###########################################################################

  //Dipole_Parameter::Reset();

#ifdef PARAMINIT_OUTPUT
  Dipole_Parameter::Show();
#endif

  extern Run_Parameter ATOOLS::rpa;
  Dipole_Parameter::Check(rpa.gen.Ecms());
#ifdef PARAMINIT_OUTPUT
  cout<<"\n  Ecms="<<rpa.gen.Ecms()<<"\n";
#endif

  Sudakov_Calculator::AdjustParameters();
  Dipole_Handler::AdjustCalcBox();
  Chain_Handler::AdjustParameters();

#ifdef PARAMINIT_OUTPUT
  Sudakov_Calculator::ShowParameters();
  Dipole_Handler::ShowCalcBox();
  Chain_Handler::ShowParameters();
#endif

  return true;

}



//=============================================================================





//eof
