#include "Basic_Func.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace std;

Kabbala Basic_Yfunc::Y(const int& z) 
{
  if (AMATOOLS::IsZero(coupl[2*z]) && AMATOOLS::IsZero(coupl[2*z+1])) {
#ifdef Kabbala_on
    return Kabbala(string("1"),Complex(1.,0.));
#else
    return Complex(1.,0.);
#endif      
  }
  //for testing purpose
  int sarg[4];

  for (short int i=0;i<4;i++) sarg[i] = arg[4*z+i].numb;

  //sarg --> &arg[4*z]

  return sgen->Get_Ynumber(sarg,&coupl[2*z],
			   Ycalc(arg[4*z].numb,arg[4*z+1].numb,
				 arg[4*z+2].numb,arg[4*z+3].numb,
				 coupl[2*z],coupl[2*z+1]));
}

template <>
Complex  Basic_Yfunc::YT<+1,+1>(const int& t1,const int& t2,const Complex& cR,const Complex& cL)
{ return BS->mu(t1)*cR*BS->eta(t2)+BS->mu(t2)*cL*BS->eta(t1); }

template <>
Complex  Basic_Yfunc::YT<-1,-1>(const int& t1,const int& t2,const Complex& cR,const Complex& cL)
{ return BS->mu(t1)*cL*BS->eta(t2)+BS->mu(t2)*cR*BS->eta(t1); }

template <>
Complex  Basic_Yfunc::YT<+1,-1>(const int& t1,const int& t2,const Complex& cR,const Complex& cL)
{ return cL*BS->S0(t1,t2); }

template <>
Complex  Basic_Yfunc::YT<-1,+1>(const int& t1,const int& t2,const Complex& cR,const Complex& cL)
{ return cR*BS->S1(t1,t2); }

Complex Basic_Yfunc::Ycalc(const int& t1,const int& sign1,
			   const int& t2,const int& sign2,
			   const Complex& cR,const Complex& cL)
{
  int sum = sign1+sign2;

  if (sum==2)  return YT<+1,+1>(t1,t2,cR,cL);
  if (sum==-2) return YT<-1,-1>(t1,t2,cR,cL);

  if ((sum==0) && (sign1==1)) return YT<+1,-1>(t1,t2,cR,cL);
  if ((sum==0) && (sign2==1)) return YT<-1,+1>(t1,t2,cR,cL);

  return Complex(0.,0.);
}

