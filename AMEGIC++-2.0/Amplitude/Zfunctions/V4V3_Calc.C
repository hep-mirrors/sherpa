#include "Calculator.H"
#include "String_Generator.H"


using namespace AMEGIC;
using namespace std;

V4V3_Calc::V4V3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS) 
{ 
  type     = zl::V4V3;
  ncoupl=18;narg=10;pn=6;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gluon4));
  lorentzlist.push_back(Lorentz_Function(lf::Gauge3));
  for (short int i=0;i<5;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[5].SetParticleArg(0,1,5,2);           
  lorentzlist[6].SetParticleArg(-5,3,4);     
}

Kabbala V4V3_Calc::ZXX(const int& a,const int& b,
		       const int& c,const int& d)
{
  return 
    -Z(b,a)*(X(c,b)+X(c,d))
    +Z(c,a)*(X(b,c)+X(b,d))
    +Z(c,b)*(X(a,b)-X(a,c));
}    

Kabbala V4V3_Calc::G4G3()
{ 
  return Z(2,0)*ZXX(1,3,4,5)+Z(1,0)*ZXX(2,4,3,5);
}

Kabbala V4V3_Calc::Do() 
{
  Kabbala factor = sgen->Get_Enumber(coupl[10])*sgen->Get_Enumber(coupl[11]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3)) &&
      IsZero(M(4)) &&
      IsZero(M(5))){

    return -factor*G4G3();
  }
  
  cerr<<"Massive V4V3"<<endl;
  abort();
  return Kabbala(string("0"),Complex(0.,0.));
}



















