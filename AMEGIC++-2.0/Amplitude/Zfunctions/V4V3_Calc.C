#include "Zfunc_Calc.H"
#include "String_Generator.H"


using namespace AMEGIC;
using namespace std;

Kabbala V4V3_Calc::ZXX(const int& a,const int& b,
		       const int& c,const int& d)
{
  return 
    -Z(b,a)*(X(c,b)+X(c,d))
    +Z(c,a)*(X(b,c)+X(b,d))
    +Z(c,b)*(X(a,b)-X(a,c));
}    

Kabbala V4V3_Calc::G4G3()
{ 
  return Z(2,0)*ZXX(1,3,4,5)+Z(1,0)*ZXX(2,4,3,5);
}

Kabbala V4V3_Calc::Do() 
{
  Kabbala factor = sgen->Get_Enumber(coupl[10])*sgen->Get_Enumber(coupl[11]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3)) &&
      IsZero(M(4)) &&
      IsZero(M(5))){

    return -factor*G4G3();
  }
  
  cerr<<"Massive V4V3"<<endl;
  abort();
  return Kabbala(string("0"),Complex(0.,0.));
}



















