#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace std;

Kabbala V5_Calc::G5()
{ 

  return
    +ZXXX1(4,1,3,2,0,5,6)  
    -ZXXX1(4,2,3,1,0,5,6)  
    +ZXXX1(3,2,4,1,0,5,6)  
    -ZXXX1(3,1,4,2,0,5,6)  

    +ZXXX2(2,1,5,3,0,4,6)
    +ZXXX2(1,2,5,4,0,3,6)
    +ZXXX2(4,3,6,2,0,1,5)
    +ZXXX2(3,4,6,1,0,2,5)
    +ZXXX3(3,4,6,5,1,2,0)
    -ZXXX3(4,3,6,5,1,2,0)
    +ZXXX3(2,1,5,6,3,4,0)
    -ZXXX3(1,2,5,6,3,4,0)

    //============================================================
    +Z(4,3)*(V(0,3)-V(0,4)-V(3,5)+V(4,5))*(-Z(1,0)*(X(2,1)+X(2,5))
					   +Z(2,0)*(X(1,2)+X(1,5)))   

    +Z(2,1)*(V(0,1)-V(0,2)-V(1,6)+V(2,6))*(-Z(4,0)*(X(3,4)+X(3,6))
					   +Z(3,0)*(X(4,3)+X(4,6)))
					   
    +Z(4,3)*Z(2,1)*(+(X(0,1)-X(0,2))*(V(0,3)-V(0,4)-V(3,5)+V(4,5))
		    +(X(0,3)-X(0,4))*(-V(0,1)+V(0,2)+V(1,6)-V(2,6))
		    +(X(0,5)-X(0,6))*(V(1,3)-V(1,4)-V(2,3)+V(2,4)));  
}

Kabbala V5_Calc::Do() 
{
  Kabbala factor = P(pn)*sgen->Get_Enumber(coupl[10])*sgen->Get_Enumber(coupl[11])
    *sgen->Get_Enumber(coupl[12]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3)) &&
      IsZero(M(4)) &&
      IsZero(M(5)) &&
      IsZero(M(6))) {
    return factor*G5();
  }
  
  cerr<<"Massive V5"<<endl;
  abort();
  return Kabbala(string("0"),Complex(0.,0.));
}

Kabbala V5_Calc::ZXXX1(const int& a,const int& b,
		       const int& c,const int& d,
		       const int& e,const int& f,
		       const int& g
		       )
{
  return Z(c,d)*(X(a,c)+X(a,g))*(X(b,d)+X(b,f))*(X(e,g)-X(e,f));
}

Kabbala V5_Calc::ZXXX2(const int& a,const int& b,
		       const int& c,const int& d,
		       const int& e,const int& f,
		       const int& g)
{
  return Z(b,e)*(X(a,b)+X(a,c))*(+(X(d,e)-X(d,c))*(X(f,d)+X(f,g))
				 -(X(f,e)-X(f,c))*(X(d,f)+X(d,g)));
}

Kabbala V5_Calc::ZXXX3(const int& a,const int& b,
		       const int& c,const int& d,
		       const int& e,const int& f,
		       const int& g)
{
  return Z(f,e)*(X(a,b)+X(a,c))*((X(b,g)-X(b,d))*(X(g,e)-X(g,f))
				+(X(b,e)-X(b,f))*(X(g,d)-X(g,c)));    
}


















