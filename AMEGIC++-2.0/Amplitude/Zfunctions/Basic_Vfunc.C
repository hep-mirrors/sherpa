#include "Basic_Func.H"
#include "Basic_Sfuncs.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;

Kabbala Basic_Vfunc::V(const int a,const int b)
{
  Complex vc = Vcalc(ps[iabs(a)].numb,ps[iabs(b)].numb);
  if ( ATOOLS::IsZero(vc) ) return sgen->Get_Enumber(Complex(0.,0.));
  
  int sign = Sign(a)*Sign(b)*ps[iabs(a)].direction*ps[iabs(b)].direction; 

  if (ps[iabs(a)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(a)].numb);
  if (ps[iabs(b)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(b)].numb);

  return (sign>0) ?
    sgen->Get_Snumber(ps[iabs(a)].numb,ps[iabs(b)].numb,vc)
    :
    -sgen->Get_Snumber(ps[iabs(a)].numb,ps[iabs(b)].numb,vc);
}

Complex Basic_Vfunc::Vcalc(const int a,const int b)
{ return BS->Momentum(a)*BS->Momentum(b);}
  
Kabbala Basic_Vfunc::Vcplx(const int a,const int b,const int s)
{
  Complex vc = Vcplxcalc(a,b);
  if ( ATOOLS::IsZero(vc) ) return sgen->Get_Enumber(Complex(0.,0.));
  
  if(s==1) return (BS->IsComplex(a)||BS->IsComplex(b)) ?
	    sgen->Get_Scplxnumber(a,b,vc)
	    : sgen->Get_Snumber(a,b,vc);
  else return (BS->IsComplex(a)||BS->IsComplex(b)) ?
	    -sgen->Get_Scplxnumber(a,b,vc)
	    : -sgen->Get_Snumber(a,b,vc);

}
    
Complex Basic_Vfunc::Vcplxcalc(const int a,const int b)
{ 
  return Complex(BS->Momentum(a)*BS->Momentum(b)-
		 BS->MomentumImg(a)*BS->MomentumImg(b),
		 BS->Momentum(a)*BS->MomentumImg(b)+
		 BS->MomentumImg(a)*BS->Momentum(b));
}




