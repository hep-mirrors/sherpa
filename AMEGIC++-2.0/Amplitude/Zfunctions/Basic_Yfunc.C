#include "Basic_Func.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace std;

Kabbala Basic_Yfunc::Y(const int z) 
{
  int sarg[4];

  for (short int i=0;i<4;i++) sarg[i] = arg[4*z+i];

  return sgen->Get_Ynumber(sarg,&coupl[2*z],
			   Ycalc(arg[4*z],arg[4*z+1],
				 arg[4*z+2],arg[4*z+3],
				 coupl[2*z],coupl[2*z+1]));
}

template <>
Complex  Basic_Yfunc::YT<+1,+1>(const int t1,const int t2,const Complex& cR,const Complex& cL)
{ return BS->Mu(t1)*cR*BS->Eta(t2)+BS->Mu(t2)*cL*BS->Eta(t1); }

template <>
Complex  Basic_Yfunc::YT<-1,-1>(const int t1,const int t2,const Complex& cR,const Complex& cL)
{ return BS->Mu(t1)*cL*BS->Eta(t2)+BS->Mu(t2)*cR*BS->Eta(t1); }

template <>
Complex  Basic_Yfunc::YT<+1,-1>(const int t1,const int t2,const Complex& cR,const Complex& cL)
{ return cL*BS->S0d(t1,t2); }

template <>
Complex  Basic_Yfunc::YT<-1,+1>(const int t1,const int t2,const Complex& cR,const Complex& cL)
{ return cR*BS->S1d(t1,t2); }

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

