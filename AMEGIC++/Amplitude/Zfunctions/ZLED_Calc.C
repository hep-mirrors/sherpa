#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"
#include "AMEGIC++/Main/Point.H"
#include "AMEGIC++/Amplitude/Zfunc.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFTENSORCALC_GETTER(FFT_Calc,"FFT","fft calculator")

FFT_Calc::FFT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Yfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="FFT";
  ncoupl=4;narg=4;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("FFT",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
}

DEFINE_ZFTENSORCALC_GETTER(VVT_Calc,"VVT","vvt calculator")

VVT_Calc::VVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVT";
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVT",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1,2);
}

DEFINE_ZFTENSORCALC_GETTER(SST_Calc,"SST","sst calculator")

SST_Calc::SST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="SST";
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("SST",LF_Key()));
  lorentzlist[0]->SetParticleArg(0,1,2);
}

DEFINE_ZFTENSORCALC_GETTER(FFVT_Calc,"FFVT","ffvt calculator")

FFVT_Calc::FFVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Zfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="FFVT";
  ncoupl=7;narg=6;pn=2;
  lorentzlist.push_back(LF_Getter::GetObject("FFVT",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist[0]->SetParticleArg(0,1);
  lorentzlist[1]->SetParticleArg(0);
}

DEFINE_ZFTENSORCALC_GETTER(VVVT_Calc,"VVVT","vvvt calculator")

VVVT_Calc::VVVT_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVVT";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVVT",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2,3);     
}

DEFINE_ZFTENSORCALC_GETTER(SSST_Calc,"SSST","ssst calculator")

SSST_Calc::SSST_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="SSST";
  ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("SSST",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);     
}

DEFINE_ZFTENSORCALC_GETTER(FFGS_Calc,"FFGS","ffgs calculator")

FFGS_Calc::FFGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Yfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS) 
{ 
  type="FFGS";
  ncoupl=3;narg=2;pn=3;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("FFGS",LF_Key()));
}

DEFINE_ZFTENSORCALC_GETTER(VVGS_Calc,"VVGS","vvgs calculator")

VVGS_Calc::VVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVGS";
  ncoupl=8;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVGS",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1,2);
}

DEFINE_ZFTENSORCALC_GETTER(SSGS_Calc,"SSGS","ssgs calculator")

SSGS_Calc::SSGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="SSGS";
  ncoupl=6;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("SSGS",LF_Key()));
  lorentzlist[0]->SetParticleArg(0,1);
}

DEFINE_ZFTENSORCALC_GETTER(FFVGS_Calc,"FFVGS","ffvgs calculator")

FFVGS_Calc::FFVGS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS),
  Basic_Mfunc(_sgen,_BS),
  Basic_Vfunc(_sgen,_BS) { 
  type="FFVGS";
  ncoupl=4;narg=4;pn=1;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFVGS",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(0);
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

void FFT_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		       Point *const p,Point *const pf,Point *&pb,
		       int *lfnumb,int *canumb)
{
  zfc->Set_Tensor(zf,p);
  if (pf==0) zfc->Set_Out(zf,0,pb,p);
  //else zfc->Set_In(zf,0,p,pf,pb);
  zfc->Set_FermionProp(zf,p,pf);
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

void VVT_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		       Point *const p,Point *const pf,Point *&pb,
		       int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zfc->Set_Tensor(zf,p);

  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
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

void SST_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		       Point *const p,Point *const pf,Point *&pb,
		       int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zfc->Set_Tensor(zf,p);
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
  zfc->SetScalarArgs(zf,scnt,pb->left);
  zfc->SetScalarArgs(zf,scnt,pb->right);
  zfc->SetScalarArgs(zf,scnt,pb->middle);
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

void FFVT_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  zfc->Set_Tensor(zf,p);
  if(pf==0){
    zfc->Set_Out(zf,0,pb,p);
    if(pb->fl.IsVector()) zfc->Set_In(zf,1,p,0,pb);
    else if(pb->left->fl.IsVector()) zfc->Set_Out(zf,1,pb->left,p);
    else if(pb->right->fl.IsVector()) zfc->Set_Out(zf,1,pb->right,p);
    else if(pb->middle->fl.IsVector()) zfc->Set_Out(zf,1,pb->middle,p);
  }
  else {
    if(!p->middle && pb->fl.IsVector()){
      zfc->Set_Out(zf,0,pb,p);
      zfc->Set_In(zf,1,p,pf,pb);
    }
    else zfc->Set_Out(zf,1,pb,p);
  }
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

void VVVT_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zfc->Set_Tensor(zf,p);

  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}

Kabbala SSST_Calc::Do() 
{
  int sarg[2];
  sarg[0]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[1]=BS->GetPolNumber(arg[2],arg[3],GetPMass(arg[2],arg[3]));
  return sgen->GetEnumber(coupl[0])* Vcplx(sarg[0],sarg[1]);     
}

void SSST_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  zfc->Set_Tensor(zf,p);

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
  zfc->SetScalarArgs(zf,scnt,pb->left);
  zfc->SetScalarArgs(zf,scnt,pb->right);
  zfc->SetScalarArgs(zf,scnt,pb->middle);
}

Kabbala FFGS_Calc::Do() 
{
  return X(0,1)-X(0,2)-Y(0) * sgen->GetEnumber(coupl[2]);
}

void FFGS_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  zfc->Set_FermionProp(zf,p,pf);
  zf->p_couplings[2]=p->cpl[2];
  if (pf==0) zfc->Set_Out(zf,0,pb,p);
  else       zfc->Set_In(zf,0,p,pf,pb);

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
}

Kabbala VVGS_Calc::Do() 
{
  return sgen->GetEnumber(coupl[6])*
    ( sgen->GetEnumber(coupl[7])*Z(1,0) +
      X(0,0)*X(1,2) + X(0,2)*X(1,1) );
}

void VVGS_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
}

Kabbala SSGS_Calc::Do() 
{
  return sgen->GetEnumber(coupl[4]) *
    ( V(0,1) + sgen->GetEnumber(coupl[5]) );
}

void SSGS_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			Point *const p,Point *const pf,Point *&pb,
			int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;
  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
  zfc->SetScalarArgs(zf,scnt,pb->left);
  zfc->SetScalarArgs(zf,scnt,pb->right);
  zfc->SetScalarArgs(zf,scnt,pb->middle);
}

Kabbala FFVGS_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
}

void FFVGS_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
			 Point *const p,Point *const pf,Point *&pb,
			 int *lfnumb,int *canumb)
{
  if(pf==0){
    zfc->Set_Out(zf,0,pb,p);
    if(pb->fl.IsVector()) zfc->Set_In(zf,1,p,0,pb);
    else if(pb->left->fl.IsVector()) zfc->Set_Out(zf,1,pb->left,p);
    else if(pb->right->fl.IsVector()) zfc->Set_Out(zf,1,pb->right,p);
    else if(pb->middle->fl.IsVector()) zfc->Set_Out(zf,1,pb->middle,p);
  }
  else {
    if(!p->middle && pb->fl.IsVector()){
      zfc->Set_Out(zf,0,pb,p);
      zfc->Set_In(zf,1,p,pf,pb);
    }
    else zfc->Set_Out(zf,1,pb,p);
  }

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
  if(pb->fl.IsVector()) pb=p;   
  zfc->SetScalarArgs(zf,scnt,pb->left);
  zfc->SetScalarArgs(zf,scnt,pb->right);
  zfc->SetScalarArgs(zf,scnt,pb->middle);
}

