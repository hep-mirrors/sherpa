#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

Y_Calc::Y_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) 
  : Basic_Func(_sgen,_BS), Zfunc_Calc(_sgen,_BS), Basic_Yfunc(_sgen,_BS) { 
  type     = zl::Y;
  ncoupl=2;narg=2;pn=1;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::FFS));
}

Z_Calc::Z_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type     = zl::Z;
  ncoupl=4;narg=4;pn=1;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(0);
}

VVS_Calc::VVS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVS;
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=5;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gab));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,1);
}

VVSS4_Calc::VVSS4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVSS4;
  ncoupl=5;narg=4;pn=2;
#ifdef Scalar_Args
  narg=6;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::VVSS));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,1);
}

SSV_Calc::SSV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::SSV;
  ncoupl=7;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::SSV));
  lorentzlist[0].SetParticleArg(2);
  lorentzlist[1].SetParticleArg(0,1,2);
}

SSS_Calc::SSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS),
  Zfunc_Calc(_sgen,_BS)
{ 
  type     = zl::SSS;
  ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
  narg=3;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::SSS));
}

SSSS_Calc::SSSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS),
  Zfunc_Calc(_sgen,_BS)
{ 
  type     = zl::SSSS;
  ncoupl=1;narg=0;pn=0;
#ifdef Scalar_Args
  narg=4;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::SSSS));
}

VVSS_Calc::VVSS_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::VVSS;
  ncoupl=6;narg=4;pn=3;
#ifdef Scalar_Args
  narg=6;
#endif
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gab));
  lorentzlist.push_back(Lorentz_Function(lf::Gab));
  lorentzlist[0].SetParticleArg(0);
  lorentzlist[1].SetParticleArg(1);
  lorentzlist[2].SetParticleArg(0,2);
  lorentzlist[3].SetParticleArg(2,1);
}

Kabbala Y_Calc::Do() {return Y(0);}

Kabbala Z_Calc::Do() 
{
  if (IsZero(M(0))) return Z(0,1);
  return (Z(0,1)-M(0)*X(0,0)*X(1,0));
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

