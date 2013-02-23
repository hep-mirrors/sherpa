#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFCALC_GETTER(G4_Calc,"G4","g4 calculator")

G4_Calc::G4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="G4";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gluon4",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala G4_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8]);
 
  return factor*(Z(0,1)*Z(2,3)-Z(0,3)*Z(2,1));
  
}



















