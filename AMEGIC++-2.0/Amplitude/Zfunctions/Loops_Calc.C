#include "Basic_Sfuncs.H"
#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;

Triangle_Calc::Triangle_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::Triangle;
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Triangle));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,1);
}

Kabbala Triangle_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[4]);
  return prefactor*(X(0,1)*X(1,0)-V(0,1)*Z(1,0));
//   return prefactor*(X(0,0)*X(1,1)+X(0,1)*X(1,0)-V(0,1)*Z(1,0));
}

Box_Calc::Box_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::Box;
  ncoupl=10;narg=6;pn=3;
#ifdef Scalar_Args
  narg=7;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Box));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(2);
  lorentzlist[3].SetParticleArg(0,1,2);
}

Kabbala Box_Calc::Do() 
{
//   Kabbala prefactor = sgen->GetEnumber(0.);
  Kabbala prefactor = sgen->GetEnumber(coupl[6]);
  return prefactor*(Z(1,0)*(X(2,0)-X(2,1))+Z(2,0)*(X(1,2)-X(1,0))+Z(2,1)*(X(0,1)-X(0,2)));
}

Pentagon_Calc::Pentagon_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::Pentagon;
  ncoupl=11;narg=8;pn=5;
#ifdef Scalar_Args
  narg=9;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gluon4));
  lorentzlist.push_back(Lorentz_Function(lf::C4GS));
  for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[4].SetParticleArg(0,1,2,4);     
  lorentzlist[5].SetParticleArg(-4,3);     
}

Kabbala Pentagon_Calc::Do() 
{
  Kabbala prefactor = sgen->GetEnumber(coupl[9]);
  return prefactor*(Z(0,1)*Z(2,3)-Z(0,3)*Z(2,1));
}

