#include "Zfunc_Calc.H"
#include "String_Generator.H"

using namespace AMEGIC;

Kabbala FFT_Calc::Do() 
{
  int sarg[4];
  //cout<<"FFT_Calc:"<<arg[4].numb<<","<<arg[5].numb<<","<<arg[6].numb<<","<<arg[7].numb<<endl;
  sarg[2]=BS->Get_Pol_Number(arg[4].numb,arg[5].numb,GetPMass(arg[4].numb,arg[5].numb));
  sarg[3]=BS->Get_Pol_Number(arg[6].numb,arg[7].numb,GetPMass(arg[6].numb,arg[7].numb));

  int s1=ps[1].direction,s2=ps[2].direction;
 
  if (ps[1].numb<BS->GetNmomenta()) s1 *= BS->Sign(ps[1].numb);
  if (ps[2].numb<BS->GetNmomenta()) s2 *= BS->Sign(ps[2].numb);
  //cout<<"FFT_Calc:"<<ps[1].numb<<","<<ps[2].numb<<":"<<sarg[2]<<","<<sarg[3]<<":"<<s1<<","<<s2<<endl;
  return 
    ( X(0,1,0)* ( Vcplx(ps[1].numb,sarg[3],s1)-Vcplx(ps[2].numb,sarg[3],s2) ) +
      X(0,1,1)* ( Vcplx(ps[1].numb,sarg[2],s1)-Vcplx(ps[2].numb,sarg[2],s2) ) -
      sgen->Get_Enumber(2.) * Vcplx(sarg[2],sarg[3]) *
      (X(0,1)-X(0,2)-Y(0) * sgen->Get_Enumber(coupl[2])) );
}

Kabbala VVT_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=BS->Get_Pol_Number(arg[8].numb,arg[9].numb,GetPMass(arg[8].numb,arg[9].numb));
  sarg[3]=BS->Get_Pol_Number(arg[10].numb,arg[11].numb,GetPMass(arg[10].numb,arg[11].numb));
  int s0=ps[0].direction,
      s1=ps[1].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
  //cout<<"VVT_Calc: "<<ps[0].numb<<":"<<ps[0].direction<<";"<<ps[1].numb<<":"<<ps[1].direction<<endl;
  Kabbala Z01 = Z(0,1),
          X01 = X(0,1)+X(0,0),
          X10 = X(1,0)+X(1,1),
          X02 = X(0,2,0),
          X03 = X(0,2,1),
          X12 = X(1,2,0),
          X13 = X(1,2,1),
          V02 = Vcplx(sarg[0],sarg[2],s0),
          V03 = Vcplx(sarg[0],sarg[3],s0),
          V12 = Vcplx(sarg[1],sarg[2],s1),
          V13 = Vcplx(sarg[1],sarg[3],s1),
          V23 = Vcplx(sarg[2],sarg[3]);
  return sgen->Get_Enumber(coupl[4])*
    ( ( sgen->Get_Enumber(coupl[5]) + V(0,1) )*( X02*X13 + X03*X12 - Z01*V23 )
      + X01*X10*V23 - ( X01*( X13*V02 + X12*V03 ) +
			X10*( X02*V13 + X03*V12 ) -
			Z01*( V03*V12 + V02*V13 ) ) );
}

Kabbala SST_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=BS->Get_Pol_Number(arg[8].numb,arg[9].numb,GetPMass(arg[8].numb,arg[9].numb));
  sarg[3]=BS->Get_Pol_Number(arg[10].numb,arg[11].numb,GetPMass(arg[10].numb,arg[11].numb));
  int s0=ps[0].direction,
      s1=ps[1].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
  return sgen->Get_Enumber(coupl[4])*
    ( ( sgen->Get_Enumber(coupl[5]) + V(0,1) ) * Vcplx(sarg[2],sarg[3]) -
      Vcplx(sarg[0],sarg[2],s0) * Vcplx(sarg[1],sarg[3],s1) -
      Vcplx(sarg[0],sarg[3],s0) * Vcplx(sarg[1],sarg[2],s1) );     
}

Kabbala FFVT_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[2]=BS->Get_Pol_Number(arg[8].numb,arg[9].numb,GetPMass(arg[8].numb,arg[9].numb));
  sarg[3]=BS->Get_Pol_Number(arg[10].numb,arg[11].numb,GetPMass(arg[10].numb,arg[11].numb));
  int s0=ps[0].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);

  Kabbala V23x2=sgen->Get_Enumber(2.) * Vcplx(sarg[2],sarg[3]);
  
  if(IsZero(M(0))) return  X(0,2,1)*X(1,2,0) + X(0,2,0)*X(1,2,1) - Z(0,1)*V23x2;
  return 
    ( X(0,2,1)*X(1,2,0) + X(0,2,0)*X(1,2,1) - Z(0,1)*V23x2
      - M(0)*( X(0,2,1)*Vcplx(sarg[0],sarg[2],s0) + X(0,2,0)*Vcplx(sarg[0],sarg[3],s0)
	       - X(0,0)*V23x2 )*X(1,0) );
}

Kabbala VVVT_Calc::Do() 
{
  int sarg[5];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->Get_Pol_Number(arg[12].numb,arg[13].numb,GetPMass(arg[12].numb,arg[13].numb));
  sarg[4]=BS->Get_Pol_Number(arg[14].numb,arg[15].numb,GetPMass(arg[14].numb,arg[15].numb));
  int s0=ps[0].direction,s1=ps[1].direction,s2=ps[2].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
  if (sarg[2]<BS->GetNmomenta()) s2 *= BS->Sign(sarg[2]);
  //cout<<"VVVT_Calc: "<<ps[0].numb<<":"<<ps[0].direction<<";"<<ps[1].numb<<":"<<ps[1].direction<<";"<<ps[2].numb<<":"<<ps[2].direction<<endl;

  Kabbala V34=Vcplx(sarg[3],sarg[4]);

  return sgen->Get_Enumber(coupl[6])*
    (  (X(2,0)-X(2,1))*(X(0,3,0)*X(1,3,1)+X(0,3,1)*X(1,3,0)-Z(0,1)*V34) +
       (X(0,1)-X(0,2))*(X(1,3,0)*X(2,3,1)+X(1,3,1)*X(2,3,0)-Z(1,2)*V34) +
       (X(1,2)-X(1,0))*(X(2,3,0)*X(0,3,1)+X(2,3,1)*X(0,3,0)-Z(2,0)*V34) +
       Z(0,1)*(X(2,3,0)*(Vcplx(sarg[0],sarg[4],s0)-Vcplx(sarg[1],sarg[4],s1)) +
	       X(2,3,1)*(Vcplx(sarg[0],sarg[3],s0)-Vcplx(sarg[1],sarg[3],s1))) +
       Z(1,2)*(X(0,3,0)*(Vcplx(sarg[1],sarg[4],s1)-Vcplx(sarg[2],sarg[4],s2)) +
	       X(0,3,1)*(Vcplx(sarg[1],sarg[3],s1)-Vcplx(sarg[2],sarg[3],s2))) +
       Z(2,0)*(X(1,3,0)*(Vcplx(sarg[2],sarg[4],s2)-Vcplx(sarg[0],sarg[4],s0)) +
	       X(1,3,1)*(Vcplx(sarg[2],sarg[3],s2)-Vcplx(sarg[0],sarg[3],s0))) );
}

Kabbala SSST_Calc::Do() 
{
  int sarg[2];
  sarg[0]=BS->Get_Pol_Number(arg[0].numb,arg[1].numb,GetPMass(arg[0].numb,arg[1].numb));
  sarg[1]=BS->Get_Pol_Number(arg[2].numb,arg[3].numb,GetPMass(arg[2].numb,arg[3].numb));
  return sgen->Get_Enumber(coupl[0])* Vcplx(sarg[0],sarg[1]);     
}

Kabbala FFGS_Calc::Do() 
{
  return X(0,1)-X(0,2)-Y(0) * sgen->Get_Enumber(coupl[2]);
}

Kabbala VVGS_Calc::Do() 
{
  return sgen->Get_Enumber(coupl[6])*
    ( sgen->Get_Enumber(coupl[7])*Z(1,0) +
      X(0,0)*X(1,2) + X(0,2)*X(1,1) );
}

Kabbala SSGS_Calc::Do() 
{
  return sgen->Get_Enumber(coupl[4]) *
    ( V(0,1) + sgen->Get_Enumber(coupl[5]) );
}

Kabbala FFVGS_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
}


