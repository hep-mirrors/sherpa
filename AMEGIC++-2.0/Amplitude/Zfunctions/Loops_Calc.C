#include "Basic_Sfuncs.H"
#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFLOOPCALC_GETTER(Triangle_Calc,TriangleCalc_Getter,"Triangle","triangle calculator")

Triangle_Calc::Triangle_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="Triangle";
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Triangle",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1);
}

Kabbala Triangle_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[4]);
  return prefactor*(X(0,1)*X(1,0)-V(0,1)*Z(1,0));
//   return prefactor*(X(0,0)*X(1,1)+X(0,1)*X(1,0)-V(0,1)*Z(1,0));
}

DEFINE_ZFLOOPCALC_GETTER(Box_Calc,BoxCalc_Getter,"Box","box calculator")

Box_Calc::Box_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="Box";
  ncoupl=10;narg=6;pn=3;
#ifdef Scalar_Args
  narg=7;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Box",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(2);
  lorentzlist[3]->SetParticleArg(0,1,2);
}

Kabbala Box_Calc::Do() 
{
//   Kabbala prefactor = sgen->GetEnumber(0.);
  Kabbala prefactor = sgen->GetEnumber(coupl[6]);
  return prefactor*(Z(1,0)*(X(2,0)-X(2,1))+Z(2,0)*(X(1,2)-X(1,0))+Z(2,1)*(X(0,1)-X(0,2)));
}

DEFINE_ZFLOOPCALC_GETTER(Pentagon_Calc,PentagonCalc_Getter,"Pentagon","pentagon calculator")

Pentagon_Calc::Pentagon_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="Pentagon";
  ncoupl=11;narg=8;pn=5;
#ifdef Scalar_Args
  narg=9;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gluon4",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("C4GS",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,4);     
  lorentzlist[5]->SetParticleArg(-4,3);     
}

Kabbala Pentagon_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[9]);
  return prefactor*(Z(0,1)*Z(2,3)-Z(0,3)*Z(2,1));
}

