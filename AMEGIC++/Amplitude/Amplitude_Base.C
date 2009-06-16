#include "AMEGIC++/Amplitude/Amplitude_Base.H"

using namespace AMEGIC;

void Amplitude_Base::SetStringOn()  
{ buildstring = 1; }

void Amplitude_Base::SetStringOff()  
{ buildstring = 0; }

Point*  Amplitude_Base::GetPointlist() 
{std::cerr<<"Error: Virtual  Amplitude_Base::GetPointlist() called!"<<std::endl;return 0;}

void Amplitude_Base::Add(Amplitude_Base* a, int sign)
{std::cerr<<"Error: Virtual  Amplitude_Base::Add() called!"<<std::endl;}

int Amplitude_Base::Size()
{
  return 1;
}
 
Zfunc_List* Amplitude_Base::GetZlist(){
  std::cerr<<"Error: Virtual  Amplitude_Base::GetZlist() called!"<<std::endl;return 0;
}

Pfunc_List* Amplitude_Base::GetPlist(){
  std::cerr<<"Error: Virtual  Amplitude_Base::GetPlist() called!"<<std::endl;return 0;
} 

int Amplitude_Base::GetSign() {
  std::cerr<<"Error: Virtual  Amplitude_Base::GetSign() called!"<<std::endl;return 0;
}

void Amplitude_Base::SetSign(int) {
  std::cerr<<"Error: Virtual  Amplitude_Base::SetSign() called!"<<std::endl;
}

void Amplitude_Base::BuildGlobalString(int* i,int j,Basic_Sfuncs* BS ,ATOOLS::Flavour* fl,
				       String_Handler* shand)
{std::cerr<<"Error: Virtual  Amplitude_Base::BuildGlobalString() called!"<<std::endl;}
