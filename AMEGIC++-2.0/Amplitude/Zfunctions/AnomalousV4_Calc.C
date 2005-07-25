#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;

AnomalousV4_Calc::AnomalousV4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::AV4;
  ncoupl=10;narg=8;pn=4;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::AGauge4));
  for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[4].SetParticleArg(0,1,2,3);     
}

Kabbala AnomalousV4_Calc::Do() 
{
  Kabbala factor1 = sgen->GetEnumber(coupl[8]);
  Kabbala factor2 = sgen->GetEnumber(coupl[9]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3))) {
    return 2*factor1*Z(0,1)*Z(2,3)-factor2*(Z(0,2)*Z(1,3)+Z(0,3)*Z(1,2));
  }
  if (!IsZero(X(0,0)*M(0)) || 
      !IsZero(X(1,1)*M(1)) || 
      !IsZero(X(2,2)*M(2)) || 
      !IsZero(X(3,3)*M(3))) {
    std::cerr<<"Error in AnomalousV4_Calc::Do(): not cutted massive vertex!"<<std::endl;
    abort();
  }
  return 2*factor1*Z(0,1)*Z(2,3)-factor2*(Z(0,2)*Z(1,3)+Z(0,3)*Z(1,2));

}














