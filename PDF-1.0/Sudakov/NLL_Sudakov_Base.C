#include "NLL_Sudakov_Base.H"

using namespace SHERPA;

NLL_Sudakov_Base::~NLL_Sudakov_Base() 
{
}

void NLL_Sudakov_Base::AssignKey(ATOOLS::Integration_Info *const info)
{
}

double NLL_Dummy_Sudakov::operator()(double,double)                 
{
  return 1.; 
}

double NLL_Dummy_Sudakov::operator()(double)                        
{ 
  return 1.; 
}

double NLL_Dummy_Sudakov::Log(double,double) 
{ 
  return 0; 
}


