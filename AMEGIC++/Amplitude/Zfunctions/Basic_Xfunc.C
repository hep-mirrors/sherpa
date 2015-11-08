#include "AMEGIC++/Amplitude/Zfunctions/Basic_Func.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/String/String_Generator.H"
#include "ATOOLS/Org/Message.H"

using namespace AMEGIC;
using namespace ATOOLS;

Kabbala Basic_Xfunc::X(const int a, const int b)
{
  int sarg[5];

  for (short int i=0;i<2;i++) sarg[i] = arg[4*a+i];
  sarg[2] = ps[iabs(b)].numb;
  for (short int i=3;i<5;i++) sarg[i] = arg[4*a+i-1];
  
  int sign = Sign(b)*ps[iabs(b)].direction;
  if (ps[iabs(b)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(b)].numb);

  //Marker for -99
  if(sarg[0]==99){
    if(sarg[3]==sarg[2]){
      return sgen->GetEnumber(Complex(0.,0.));
    }
    return (sign>0) ? 
      Vcplx(BS->GetPolNumber(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]) 
      :
      -Vcplx(BS->GetPolNumber(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]); 
  }
  //Marker for -99
  if(sarg[3]==99){
    if(sarg[0]==sarg[2]){
      return sgen->GetEnumber(Complex(0.,0.));
    }
    return (sign>0) ? 
      Vcplx(BS->GetPolNumber(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2])
      :
      -Vcplx(BS->GetPolNumber(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2]);
  }

  return (sign>0) ? 
    sgen->GetXnumber(sarg,&coupl[2*a],
		      Xcalc(arg[4*a]  ,arg[4*a+1],
			    ps[iabs(b)].numb,
			    arg[4*a+2],arg[4*a+3],
			    coupl[2*a],coupl[2*a+1]))
    :
    -sgen->GetXnumber(sarg,&coupl[2*a],
		       Xcalc(arg[4*a]  ,arg[4*a+1],
			     ps[iabs(b)].numb,
			     arg[4*a+2],arg[4*a+3],
			     coupl[2*a],coupl[2*a+1]));
}

Kabbala Basic_Xfunc::X(const int a, const int b, const int m)
{
  int sarg[5];
  
  for (short int i=0;i<2;i++) sarg[i] = arg[4*a+i];
  sarg[2] = BS->GetPolNumber(arg[4*b+2*m],arg[4*b+2*m+1],
			     GetPMass(arg[4*b+2*m],arg[4*b+2*m+1]));
  for (short int i=3;i<5;i++) sarg[i] = arg[4*a+i-1];
  
  //Marker for -99
  if(sarg[0]==99){
    if(sarg[3]==sarg[2]){
      return sgen->GetEnumber(Complex(0.,0.));
    }
    return Vcplx(BS->GetPolNumber(sarg[3],sarg[4],GetPMass(sarg[3],sarg[4])),sarg[2]); 
  }
  //Marker for -99
  if(sarg[3]==99){
    if(sarg[0]==sarg[2]){
      return sgen->GetEnumber(Complex(0.,0.));
    }
    return Vcplx(BS->GetPolNumber(sarg[0],sarg[1],GetPMass(sarg[0],sarg[1])),sarg[2]);
  }

  return sgen->GetXnumber(sarg,&coupl[2*a],
		      Xcalc(arg[4*a]  ,arg[4*a+1],
			    sarg[2],
			    arg[4*a+2],arg[4*a+3],
			    coupl[2*a],coupl[2*a+1]));
}


Complex Basic_Xfunc::Xcalc(const int t1,const int sign1,const int t2,
			   const int t3,const int sign3,
			   const Complex& cR,const Complex& cL)
{
  
  int sum = sign1+sign3;

  if (sum==2){
    if (BS->IsMomSum(t2,t1,t3)) 
      if(BS->Mu(t1)==Complex(0.0,0.0) && BS->Mu(t3)==Complex(0.0,0.0)) return Complex(0.,0.);
    return cL*BS->Eta(t2)*BS->Eta(t2)*BS->Mu(t1)*BS->Mu(t3)+
           cR*(BS->Eta(t1)*BS->Eta(t3)*BS->Mu(t2)*BS->Mu(t2)+
	       BS->S0(t1,t2)*BS->S1(t2,t3));
  }  
  if (sum==-2){
    if (BS->IsMomSum(t2,t1,t3)) 
      if(BS->Mu(t1)==Complex(0.0,0.0) && BS->Mu(t3)==Complex(0.0,0.0)) return Complex(0.,0.);
    return cR*BS->Eta(t2)*BS->Eta(t2)*BS->Mu(t1)*BS->Mu(t3)+
           cL*(BS->Eta(t1)*BS->Eta(t3)*BS->Mu(t2)*BS->Mu(t2)+
	       BS->S1(t1,t2)*BS->S0(t2,t3));
    
  }

  if(t1==t3) if(cL==cR) return Complex(0.,0.);

  if (sign1==1) 
    return BS->Eta(t2)*(cL*BS->Mu(t1)*BS->S0(t2,t3)+
			cR*BS->Mu(t3)*BS->S0(t1,t2));
  
  if (sign3==1) 
    return BS->Eta(t2)*(cR*BS->Mu(t1)*BS->S1(t2,t3)+
			cL*BS->Mu(t3)*BS->S1(t1,t2));

  return 0.;
}









