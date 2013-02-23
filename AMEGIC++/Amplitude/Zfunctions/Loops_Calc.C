#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"
#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFLOOPCALC_GETTER(Triangle_Calc,"Triangle","triangle calculator")

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

DEFINE_ZFLOOPCALC_GETTER(Box_Calc,"Box","box calculator")

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

DEFINE_ZFLOOPCALC_GETTER(Pentagon_Calc,"Pentagon","pentagon calculator")

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


DEFINE_ZFLOOPCALC_GETTER(PseudoTriangle_Calc,"PseudoTriangle","pseudo triangle calculator")

PseudoTriangle_Calc::PseudoTriangle_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS) 
{ 
  type="PseudoTriangle";
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("PseudoTriangle",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(0,1);
}

Kabbala PseudoTriangle_Calc::Do() 
{
  int sarg[4];

  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  int si;
  si=(arg[0]==99)*2;
  sarg[2]=BS->GetPolNumber(arg[0+si],arg[1+si],GetPMass(arg[0+si],arg[1+si]));
  si=(arg[4]==99)*2;
  sarg[3]=BS->GetPolNumber(arg[4+si],arg[5+si],GetPMass(arg[4+si],arg[5+si]));

  Kabbala prefactor = sgen->GetEnumber(coupl[4]);
  return prefactor*Epsilon(sarg[0],sarg[2],sarg[1],sarg[3],1);
}

DEFINE_ZFLOOPCALC_GETTER(PseudoBox_Calc,"PseudoBox","pseudo box calculator")

PseudoBox_Calc::PseudoBox_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS) 
{ 
  type="PseudoBox";
  ncoupl=10;narg=6;pn=3;
#ifdef Scalar_Args
  narg=7;
#endif
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("PseudoBox",LF_Key()));
  lorentzlist[0]->SetParticleArg(0);
  lorentzlist[1]->SetParticleArg(1);
  lorentzlist[2]->SetParticleArg(2);
  lorentzlist[3]->SetParticleArg(0,1,2);
}

Kabbala PseudoBox_Calc::Do() 
{
  int sarg[6];

  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  int si;
  si=(arg[0]==99)*2;
  sarg[3]=BS->GetPolNumber(arg[0+si],arg[1+si],GetPMass(arg[0+si],arg[1+si]));
  si=(arg[4]==99)*2;
  sarg[4]=BS->GetPolNumber(arg[4+si],arg[5+si],GetPMass(arg[4+si],arg[5+si]));
  si=(arg[8]==99)*2;
  sarg[5]=BS->GetPolNumber(arg[8+si],arg[9+si],GetPMass(arg[8+si],arg[9+si]));

  Kabbala prefactor = sgen->GetEnumber(coupl[6]);
  return prefactor*(Epsilon(sarg[0],sarg[3],sarg[4],sarg[5],1)+
		    Epsilon(sarg[1],sarg[4],sarg[5],sarg[3],1)+
		    Epsilon(sarg[2],sarg[5],sarg[3],sarg[4],1));
}
