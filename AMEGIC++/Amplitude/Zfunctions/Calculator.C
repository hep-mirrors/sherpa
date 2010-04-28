#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"
#include "AMEGIC++/Amplitude/Zfunc_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFCALC_GETTER(Y_Calc,YCalc_Getter,"Y","y calculator")

Y_Calc::Y_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) 
  : Basic_Func(_sgen,_BS), Zfunc_Calc(_sgen,_BS), Basic_Yfunc(_sgen,_BS) { 
  type="Y";
  ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("FFS",LF_Key()));
}

DEFINE_ZFCALC_GETTER(Z_Calc,ZCalc_Getter,"Z","z calculator")

Z_Calc::Z_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="Z";
  ncoupl=4;narg=4;pn=1;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(0);
}

DEFINE_ZFCALC_GETTER(VVS_Calc,VVSCalc_Getter,"VVS","vvs calculator")

VVS_Calc::VVS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVS";
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gab",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1);
}

DEFINE_ZFCALC_GETTER(VVSS4_Calc,VVSS4Calc_Getter,"VVSS4","vvss4 calculator")

VVSS4_Calc::VVSS4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVSS4";
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=6;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVSS",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1);
}

DEFINE_ZFCALC_GETTER(SSV_Calc,SSVCalc_Getter,"SSV","ssv calculator")

SSV_Calc::SSV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="SSV";
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("SSV",LF_Key()));
  lorentzlist[0]->SetParticleArg(2);
  lorentzlist[1]->SetParticleArg(0,1,2);
}

DEFINE_ZFCALC_GETTER(SSS_Calc,SSSCalc_Getter,"SSS","sss calculator")

SSS_Calc::SSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS),
  Zfunc_Calc(_sgen,_BS)
{ 
  type="SSS";
  ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("SSS",LF_Key()));
}

DEFINE_ZFCALC_GETTER(SSSS_Calc,SSSSCalc_Getter,"SSSS","ssss calculator")

SSSS_Calc::SSSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS),
  Zfunc_Calc(_sgen,_BS)
{ 
  type="SSSS";
  ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
  narg=4;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("SSSS",LF_Key()));
}

DEFINE_ZFCALC_GETTER(VVSS_Calc,VVSSCalc_Getter,"VVSS","vvss calculator")

VVSS_Calc::VVSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="VVSS";
  ncoupl=6;narg=4;pn=3;
#ifdef Scalar_Args
  narg=6;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gab",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gab",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,2);
  lorentzlist[3]->SetParticleArg(2,1);
}

DEFINE_ZFCALC_GETTER(SSVgen_Calc,SSVgenCalc_Getter,"SSVgen","LHTM SSV calculator")

  SSVgen_Calc::SSVgen_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) :
    Basic_Func(_sgen,_BS),
    Zfunc_Calc(_sgen,_BS),
    Basic_Xfunc(_sgen,_BS),
    Basic_Mfunc(_sgen,_BS),
    Basic_Vfunc(_sgen,_BS)
{
  type="SSVgen";
  ncoupl=8;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("SSVgen",LF_Key()));
  lorentzlist[0]->SetParticleArg(2);
  lorentzlist[1]->SetParticleArg(0,1,2);
}

Kabbala Y_Calc::Do() {return Y(0);}

void Y_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		     Point *const p,Point *const pf,Point *&pb,
		     int *lfnumb,int *canumb)
{
  if (pf==0) zfc->Set_Out(zf,0,pb,p);
  else       zfc->Set_In(zf,0,p,pf,pb);

  int scnt(narg-GetScalarNumb());
  if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
}

Kabbala Z_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
}

void Z_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
		     Point *const p,Point *const pf,Point *&pb,
		     int *lfnumb,int *canumb)
{
  zfc->Set_Out(zf,1,pb,p);
}

Kabbala VVS_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[4]);
  return prefactor*(-M(1)*X(0,1)*X(1,1)+X(0,0)*(-M(0)*X(1,0)+M(0)*M(1)*V(0,1)*X(1,1))+Z(1,0));
}

Kabbala VVSS4_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[4]);
  return prefactor*(-M(1)*X(0,1)*X(1,1)+X(0,0)*(-M(0)*X(1,0)+M(0)*M(1)*V(0,1)*X(1,1))+Z(1,0));
}

Kabbala SSV_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[6]);
  return prefactor*(-X(2,0)+X(2,1)+M(2)*(V(0,2)+V(1,2))*X(2,2));
}

Kabbala SSS_Calc::Do() {return sgen->GetEnumber(coupl[0]);}

Kabbala SSSS_Calc::Do() {return sgen->GetEnumber(coupl[0]);}

Kabbala VVSS_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[4])*sgen->GetEnumber(coupl[5]);
  return -prefactor*( M(0)*X(0,0)*X(1,0)
		     +M(1)*X(0,1)*X(1,1)
		     +M(2)*X(0,2)*X(1,2)
		     -M(0)*M(1)*V(0,1)*X(0,0)*X(1,1)
		     -M(0)*M(2)*V(0,2)*X(0,0)*X(1,2)
		     -M(1)*M(2)*V(1,2)*X(0,2)*X(1,1)
		     +M(0)*M(1)*M(2)*V(0,2)*V(1,2)*X(0,0)*X(1,1)
		     -Z(1,0));
}

Kabbala SSVgen_Calc::Do()
{
  Kabbala kcpl0 = sgen->GetEnumber(coupl[6]);
  Kabbala kcpl1 = sgen->GetEnumber(coupl[7]);
  Kabbala A = kcpl0 + kcpl1;
  Kabbala B = kcpl1 - kcpl0;

  return ( -A*X(2,0) - B*X(2,1) + M(2)*( A*V(0,2) + B*V(1,2) )*X(2,2) );
}

void SSVgen_Calc::SetArgs(Zfunc_Generator *const zfc,Zfunc *const zf,
                          Point *const p,Point *const pf,Point *&pb,
                          int *lfnumb,int *canumb)
{
  int icoupl(zf->m_narg-GetScalarNumb());

  pb->cpl.resize(2);
  zf->p_couplings[icoupl] = pb->cpl[0];icoupl++;
  zf->p_couplings[icoupl] = pb->cpl[1];icoupl++;

  zfc->SetArgs(zf,lfnumb,canumb,pb->left,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->right,p,icoupl);
  zfc->SetArgs(zf,lfnumb,canumb,pb->middle,p,icoupl);
  if(GetScalarNumb()>0){
    int scnt(narg-GetScalarNumb());
    if(pb->fl.IsScalar()) zfc->SetScalarArgs(zf,scnt,pb);
    zfc->SetScalarArgs(zf,scnt,pb->left);
    zfc->SetScalarArgs(zf,scnt,pb->right);
    zfc->SetScalarArgs(zf,scnt,pb->middle);
  }
}


