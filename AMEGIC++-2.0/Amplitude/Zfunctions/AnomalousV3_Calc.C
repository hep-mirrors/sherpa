#include "Calculator.H"
#include "String_Generator.H"
#include "Message.H"

using namespace AMEGIC;

AnomalousV3_Calc::AnomalousV3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type     = zl::AV3;
  ncoupl=10;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::AGauge3));
  for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[3].SetParticleArg(0,1,2);     
}

Kabbala AnomalousV3_Calc::Do() 
{
  if (!IsZero(X(0,0)*M(0)) || 
      !IsZero(X(1,1)*M(1)) || 
      !IsZero(X(2,2)*M(2))) {
    std::cerr<<"Error in AnomalousV3_Calc::Do(): not cutted massive vertex!"<<std::endl;
    abort();
  }

  Kabbala f1 = -sgen->GetEnumber(coupl[6]);
  Kabbala f2 = -sgen->GetEnumber(coupl[7]);
  Kabbala f3 = -sgen->GetEnumber(coupl[8]);
  Kabbala f4 = -sgen->GetEnumber(coupl[9]);

  Kabbala t1 = Z(2,1)*(X(0,1)-X(0,2));
  Kabbala t2 = (X(0,1)-X(0,2))*X(1,0)*X(2,0);
  Kabbala t3 = Z(1,0)*X(2,0)-Z(2,0)*X(1,0);
  Kabbala t4 = -Z(1,0)*X(2,0)-Z(2,0)*X(1,0);
  return f1*t1+f2*(0.5*V(0,0)*t1-t2)+f3*t3+f4*t4;
}
