#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

Kabbala Basic_Yfunc::Y(const int z) 
{
  int sarg[4];

  for (short int i=0;i<4;i++) sarg[i] = arg[4*z+i];

  return sgen->GetYnumber(sarg,&coupl[2*z],
			   Ycalc(arg[4*z],arg[4*z+1],
				 arg[4*z+2],arg[4*z+3],
				 coupl[2*z],coupl[2*z+1]));
}

Complex Basic_Yfunc::Ycalc(const int t1,const int sign1,
			   const int t2,const int sign2,
			   const Complex& cR,const Complex& cL)
{
  int sum = sign1+sign2;

  if (sum==2)  return YT<+1,+1>(t1,t2,cR,cL);
  if (sum==-2) return YT<-1,-1>(t1,t2,cR,cL);

  if ((sum==0) && (sign1==1)) return cL*BS->S0(t1,t2);
  if ((sum==0) && (sign2==1)) return cR*BS->S1(t1,t2);
  //if ((sum==0) && (sign1==1)) return YT<+1,-1>(t1,t2,cR,cL);
  //if ((sum==0) && (sign2==1)) return YT<-1,+1>(t1,t2,cR,cL);

  return Complex(0.,0.);
}

