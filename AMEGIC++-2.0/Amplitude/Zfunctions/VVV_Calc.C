#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;

Kabbala VVV_Calc::GGG() 
{ return Z(1,0)*(X(2,0)-X(2,1))+Z(2,0)*(X(1,2)-X(1,0))+Z(2,1)*(X(0,1)-X(0,2));}

Kabbala VVV_Calc::Do() 
{
  Kabbala factor = sgen->Get_Enumber(coupl[6]);

  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)))
    return factor*GGG();

  return factor*( M(0)*M(1)*V(1,2)*X(0,0)*X(1,1)*X(2,0)-M(1)*X(0,1)*X(1,1)*X(2,0)
		 -M(0)*X(0,0)*X(1,2)*X(2,0)+M(0)*X(0,0)*X(1,0)*X(2,1)
		 -M(0)*M(1)*V(0,2)*X(0,0)*X(1,1)*X(2,1)+M(1)*X(0,2)*X(1,1)*X(2,1)
		 -M(0)*M(2)*V(1,2)*X(0,0)*X(1,0)*X(2,2)+M(2)*X(0,2)*X(1,0)*X(2,2)
		 +M(1)*M(2)*V(0,2)*X(0,1)*X(1,1)*X(2,2)
		 -M(1)*M(2)*V(0,1)*X(0,2)*X(1,1)*X(2,2)
		 +M(0)*M(2)*V(0,1)*X(0,0)*X(1,2)*X(2,2)-M(2)*X(0,1)*X(1,2)*X(2,2)
		 -M(2)*V(0,2)*X(2,2)*Z(1,0)+M(2)*V(1,2)*X(2,2)*Z(1,0)
		 +M(1)*V(0,1)*X(1,1)*Z(2,0)-M(1)*V(1,2)*X(1,1)*Z(2,0)
		 -M(0)*V(0,1)*X(0,0)*Z(2,1)+M(0)*V(0,2)*X(0,0)*Z(2,1)
		 -X(1,0)*Z(2,0)+X(1,2)*Z(2,0)
		 +X(2,0)*Z(1,0)-X(2,1)*Z(1,0)
		 +X(0,1)*Z(2,1)-X(0,2)*Z(2,1));
}
