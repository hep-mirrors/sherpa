#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;

Kabbala G4_Calc::Do() 
{
  Kabbala factor = sgen->Get_Enumber(coupl[8]);
 
  return factor*(Z(0,1)*Z(2,3)-Z(0,3)*Z(2,1));
  
}



















