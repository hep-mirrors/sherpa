#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFCALC_GETTER(G4A_Calc,"VVVVA","g4 calculator")

G4A_Calc::G4A_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="VVVVA";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVVVA",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala G4A_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8]);
 
  return factor*(Z(0,3)*Z(2,1)-Z(0,2)*Z(3,1));
  
}

DEFINE_ZFCALC_GETTER(G4B_Calc,"VVVVB","g4 calculator")

G4B_Calc::G4B_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="VVVVB";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVVVB",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala G4B_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8]);
 
  return factor*(Z(0,3)*Z(1,2)-Z(0,1)*Z(3,2));
  
}

DEFINE_ZFCALC_GETTER(G4C_Calc,"VVVVC","g4 calculator")

G4C_Calc::G4C_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS),
  Basic_Xfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="VVVVC";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("FFV",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("VVVVC",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala G4C_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8]);
 
  return factor*(Z(0,2)*Z(1,3)-Z(0,1)*Z(2,3));
  
}



















