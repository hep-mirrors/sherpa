#include "Basic_Func.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace AMATOOLS;

Kabbala Basic_Vfunc::V(const int& a,const int &b)
{
  int sign = Sign(a)*Sign(b)*ps[iabs(a)].direction*ps[iabs(b)].direction; 

  if (ps[iabs(a)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(a)].numb);
  if (ps[iabs(b)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(b)].numb);

  return (sign>0) ?
    sgen->Get_Snumber(ps[iabs(a)].numb,ps[iabs(b)].numb,Vcalc(ps[iabs(a)].numb,ps[iabs(b)].numb))
    :
    -sgen->Get_Snumber(ps[iabs(a)].numb,ps[iabs(b)].numb,Vcalc(ps[iabs(a)].numb,ps[iabs(b)].numb));
}

Complex Basic_Vfunc::Vcalc(const int& a,const int &b)
{ return BS->Momentum(a)*BS->Momentum(b);}
  
Kabbala Basic_Vfunc::Vcplx(const int& a,const int &b,int s)
{
  int sarg[2];
  sarg[0]=a;
  sarg[1]=b;
  
  if(s==1) return (BS->iscplx(sarg[0])||BS->iscplx(sarg[1])) ?
	    sgen->Get_Scplxnumber(sarg[0],sarg[1],Vcplxcalc(sarg[0],sarg[1]))
	    : sgen->Get_Snumber(sarg[0],sarg[1],Vcalc(sarg[0],sarg[1]));
  else return (BS->iscplx(sarg[0])||BS->iscplx(sarg[1])) ?
	    -sgen->Get_Scplxnumber(sarg[0],sarg[1],Vcplxcalc(sarg[0],sarg[1]))
	    : -sgen->Get_Snumber(sarg[0],sarg[1],Vcalc(sarg[0],sarg[1]));

}
    
Complex Basic_Vfunc::Vcplxcalc(const int& a,const int &b)
{ 
  return Complex(BS->Momentum(a)*BS->Momentum(b)-
		 BS->Momentum_img(a)*BS->Momentum_img(b),
		 BS->Momentum(a)*BS->Momentum_img(b)+
		 BS->Momentum_img(a)*BS->Momentum(b));
}




