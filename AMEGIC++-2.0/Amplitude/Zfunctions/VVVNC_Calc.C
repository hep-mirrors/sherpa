#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;

Kabbala VVVNC_Calc::Do() 
{
  Kabbala factor = P(pn)*sgen->Get_Enumber(coupl[6])*sin(VNC(1,2));
  factor *= exp(sgen->Get_Enumber(Complex(0.,1.))*(VNC(0)+VNC(1)+VNC(2)));

 return factor*(Z(1,0)*(X(2,0)-X(2,1))+Z(2,0)*(X(1,2)-X(1,0))+Z(2,1)*(X(0,1)-X(0,2)));
}
