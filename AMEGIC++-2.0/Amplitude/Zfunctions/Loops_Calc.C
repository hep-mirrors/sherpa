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
  return prefactor*(X(0,0)*X(1,1)+X(0,1)*X(1,0)-V(0,1)*Z(1,0));
}
