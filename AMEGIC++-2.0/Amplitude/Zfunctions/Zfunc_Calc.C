#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;

Kabbala Y_Calc::Do() {return Y(0);}

Kabbala Z_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
}

Kabbala VVS_Calc::Do() 
{
  Kabbala prefactor = sgen->Get_Enumber(coupl[4]);
  return prefactor*(-M(1)*X(0,1)*X(1,1)+X(0,0)*(-M(0)*X(1,0)+M(0)*M(1)*V(0,1)*X(1,1))+Z(1,0));
}

Kabbala VVSS4_Calc::Do() 
{
  Kabbala prefactor = sgen->Get_Enumber(coupl[4]);
  return prefactor*(-M(1)*X(0,1)*X(1,1)+X(0,0)*(-M(0)*X(1,0)+M(0)*M(1)*V(0,1)*X(1,1))+Z(1,0));
}

Kabbala SSV_Calc::Do() 
{
  Kabbala prefactor = sgen->Get_Enumber(coupl[6]);
  return prefactor*(-X(2,0)+X(2,1)+M(2)*(V(0,2)+V(1,2))*X(2,2));
}

Kabbala SSS_Calc::Do() {return sgen->Get_Enumber(coupl[0]);}

Kabbala SSSS_Calc::Do() {return sgen->Get_Enumber(coupl[0]);}

Kabbala VVSS_Calc::Do() 
{
  Kabbala prefactor = sgen->Get_Enumber(coupl[4])*sgen->Get_Enumber(coupl[5]);
  return -prefactor*( M(0)*X(0,0)*X(1,0)
		     +M(1)*X(0,1)*X(1,1)
		     +M(2)*X(0,2)*X(1,2)
		     -M(0)*M(1)*V(0,1)*X(0,0)*X(1,1)
		     -M(0)*M(2)*V(0,2)*X(0,0)*X(1,2)
		     -M(1)*M(2)*V(1,2)*X(0,2)*X(1,1)
		     +M(0)*M(1)*M(2)*V(0,2)*V(1,2)*X(0,0)*X(1,1)
		     -Z(1,0));
}

//Non-Commutative QED

Kabbala ZNC_Calc::Do() 
{
  Kabbala prefactor = exp(sgen->Get_Enumber(Complex(0.,1.))*(VNC(0)+VNC(1)));
  if (IsZero(M(0))) return prefactor*Z(0,1);
  return prefactor*(Z(0,1)-M(0)*X(0,0)*X(1,0));
}
