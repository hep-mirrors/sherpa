#include "Basic_Sfuncs.H"
#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;

FFT_Calc::FFT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Yfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type     = zl::FFT;
  ncoupl=4;narg=4;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::FFT));
  lorentzlist[0].SetParticleArg(0);
}

VVT_Calc::VVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVT;
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::VVT));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,1,2);
}

SST_Calc::SST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::SST;
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::SST));
  lorentzlist[0].SetParticleArg(0,1,2);
}

FFVT_Calc::FFVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Zfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type     = zl::FFVT;
  ncoupl=7;narg=6;pn=2;
  lorentzlist.push_back(Lorentz_Function(lf::FFVT));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist[0].SetParticleArg(0,1);
  lorentzlist[1].SetParticleArg(0);
}

VVVT_Calc::VVVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVVT;
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::VVVT));
  for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[3].SetParticleArg(0,1,2,3);     
}

SSST_Calc::SSST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::SSST;
  ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::SSST));
  lorentzlist[0].SetParticleArg(0);     
}

FFGS_Calc::FFGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Yfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS) 
{ 
  type     = zl::FFGS;
  ncoupl=3;narg=2;pn=3;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::FFGS));
}


VVGS_Calc::VVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVGS;
  ncoupl=8;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::VVGS));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,1,2);
}

SSGS_Calc::SSGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::SSGS;
  ncoupl=6;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::SSGS));
  lorentzlist[0].SetParticleArg(0,1);
}

FFVGS_Calc::FFVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS),
  Basic_Mfunc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) { 
  type     = zl::FFVGS;
  ncoupl=4;narg=4;pn=1;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::FFVGS));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(0);
}


Kabbala FFT_Calc::Do() 
{
  int sarg[4];
  sarg[2]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[3]=BS->GetPolNumber(arg[6],arg[7],GetPMass(arg[6],arg[7]));

  int s1=ps[1].direction,s2=ps[2].direction;
 
  if (ps[1].numb<BS->GetNmomenta()) s1 *= BS->Sign(ps[1].numb);
  if (ps[2].numb<BS->GetNmomenta()) s2 *= BS->Sign(ps[2].numb);
  return 
    ( X(0,1,0)* ( Vcplx(ps[1].numb,sarg[3],s1)-Vcplx(ps[2].numb,sarg[3],s2) ) +
      X(0,1,1)* ( Vcplx(ps[1].numb,sarg[2],s1)-Vcplx(ps[2].numb,sarg[2],s2) ) -
      sgen->GetEnumber(2.) * Vcplx(sarg[2],sarg[3]) *
      (X(0,1)-X(0,2)-Y(0) * sgen->GetEnumber(coupl[2])) );
}

Kabbala VVT_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));
  sarg[3]=BS->GetPolNumber(arg[10],arg[11],GetPMass(arg[10],arg[11]));
  int s0=ps[0].direction,
      s1=ps[1].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
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
  return sgen->GetEnumber(coupl[4])*
    ( ( sgen->GetEnumber(coupl[5]) + V(0,1) )*( X02*X13 + X03*X12 - Z01*V23 )
      + X01*X10*V23 - ( X01*( X13*V02 + X12*V03 ) +
			X10*( X02*V13 + X03*V12 ) -
			Z01*( V03*V12 + V02*V13 ) ) );
}

Kabbala SST_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));
  sarg[3]=BS->GetPolNumber(arg[10],arg[11],GetPMass(arg[10],arg[11]));
  int s0=ps[0].direction,
      s1=ps[1].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
  return sgen->GetEnumber(coupl[4])*
    ( ( sgen->GetEnumber(coupl[5]) + V(0,1) ) * Vcplx(sarg[2],sarg[3]) -
      Vcplx(sarg[0],sarg[2],s0) * Vcplx(sarg[1],sarg[3],s1) -
      Vcplx(sarg[0],sarg[3],s0) * Vcplx(sarg[1],sarg[2],s1) );     
}

Kabbala FFVT_Calc::Do() 
{
  int sarg[4];
  sarg[0]=ps[0].numb;
  sarg[2]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));
  sarg[3]=BS->GetPolNumber(arg[10],arg[11],GetPMass(arg[10],arg[11]));
  int s0=ps[0].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);

  Kabbala V23x2=sgen->GetEnumber(2.) * Vcplx(sarg[2],sarg[3]);
  
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
  sarg[3]=BS->GetPolNumber(arg[12],arg[13],GetPMass(arg[12],arg[13]));
  sarg[4]=BS->GetPolNumber(arg[14],arg[15],GetPMass(arg[14],arg[15]));
  int s0=ps[0].direction,s1=ps[1].direction,s2=ps[2].direction;
  if (sarg[0]<BS->GetNmomenta()) s0 *= BS->Sign(sarg[0]);
  if (sarg[1]<BS->GetNmomenta()) s1 *= BS->Sign(sarg[1]);
  if (sarg[2]<BS->GetNmomenta()) s2 *= BS->Sign(sarg[2]);

  Kabbala V34=Vcplx(sarg[3],sarg[4]);

  return sgen->GetEnumber(coupl[6])*
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
  sarg[0]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[1]=BS->GetPolNumber(arg[2],arg[3],GetPMass(arg[2],arg[3]));
  return sgen->GetEnumber(coupl[0])* Vcplx(sarg[0],sarg[1]);     
}

Kabbala FFGS_Calc::Do() 
{
  return X(0,1)-X(0,2)-Y(0) * sgen->GetEnumber(coupl[2]);
}

Kabbala VVGS_Calc::Do() 
{
  return sgen->GetEnumber(coupl[6])*
    ( sgen->GetEnumber(coupl[7])*Z(1,0) +
      X(0,0)*X(1,2) + X(0,2)*X(1,1) );
}

Kabbala SSGS_Calc::Do() 
{
  return sgen->GetEnumber(coupl[4]) *
    ( V(0,1) + sgen->GetEnumber(coupl[5]) );
}

Kabbala FFVGS_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
}


