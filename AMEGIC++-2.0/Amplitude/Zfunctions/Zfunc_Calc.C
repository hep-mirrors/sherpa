#include "Zfunc_Calc.H"

using namespace AMEGIC;
using namespace ATOOLS;

#define COMPILE__Getter_Function
#define PARAMETER_TYPE AMEGIC::ZFCalc_Key
#define OBJECT_TYPE AMEGIC::Zfunc_Calc
#include "Getter_Function.C"

Zfunc_Calc::~Zfunc_Calc() 
{
  for (size_t i(0);i<lorentzlist.size();++i) delete lorentzlist[i];
}

Kabbala Zfunc_Calc::Do() 
{ 
  std::cerr<<"Error: Virtual method Zfunc_Calc::Do() called!"<<std::endl; 
  return Kabbala();
}

int Zfunc_Calc::GetScalarNumb() 
{ 
  return 0; 
}
