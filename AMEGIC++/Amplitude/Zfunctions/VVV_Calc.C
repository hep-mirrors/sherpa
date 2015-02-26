#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFCALC_GETTER(VVV_Calc,"VVV","vvv calculator")

VVV_Calc::VVV_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS)
{ 
  type="VVV";
  ncoupl=10;narg=6;pn=3;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gauge3",LF_Key()));
  for (short int i=0;i<3;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[3]->SetParticleArg(0,1,2);     
}

Kabbala VVV_Calc::GGG() 
{ return Z(1,0)*(X(2,0)-X(2,1))+Z(2,0)*(X(1,2)-X(1,0))+Z(2,1)*(X(0,1)-X(0,2));}

Kabbala VVV_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[6]);

  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)))
    return factor*GGG();

  if (!IsZero(X(0,0)*M(0)) || 
      !IsZero(X(1,1)*M(1)) || 
      !IsZero(X(2,2)*M(2))) {
    std::cerr<<"Error in VVV_Calc::Do(): not cutted massive vertex!"<<std::endl;

    return factor*( M(0)*M(1)*V(1,2)*X(0,0)*X(1,1)*X(2,0)-M(1)*X(0,1)*X(1,1)*X(2,0)
                    -M(0)*X(0,0)*X(1,2)*X(2,0)+M(0)*X(0,0)*X(1,0)*X(2,1)
                    -M(0)*M(1)*V(0,2)*X(0,0)*X(1,1)*X(2,1)+M(1)*X(0,2)*X(1,1)*X(2,1)
                    -M(0)*M(2)*V(1,2)*X(0,0)*X(1,0)*X(2,2)+M(2)*X(0,2)*X(1,0)*X(2,2)
                    +M(1)*M(2)*V(0,2)*X(0,1)*X(1,1)*X(2,2)
                    -M(1)*M(2)*V(0,1)*X(0,2)*X(1,1)*X(2,2)
                    +M(0)*M(2)*V(0,1)*X(0,0)*X(1,2)*X(2,2)-M(2)*X(0,1)*X(1,2)*X(2,2)
                    -M(2)*V(0,2)*X(2,2)*Z(1,0)+M(2)*V(1,2)*X(2,2)*Z(1,0)
                    +M(1)*V(0,1)*X(1,1)*Z(2,0)-M(1)*V(1,2)*X(1,1)*Z(2,0)
                    -M(0)*V(0,1)*X(0,0)*Z(2,1)+M(0)*V(0,2)*X(0,0)*Z(2,1)
                    -X(1,0)*Z(2,0)+X(1,2)*Z(2,0)
                    +X(2,0)*Z(1,0)-X(2,1)*Z(1,0)
                    +X(0,1)*Z(2,1)-X(0,2)*Z(2,1));
  }
  return factor*GGG();

  /*return factor*( M(0)*M(1)*V(1,2)*X(0,0)*X(1,1)*X(2,0)-M(1)*X(0,1)*X(1,1)*X(2,0)
		 -M(0)*X(0,0)*X(1,2)*X(2,0)+M(0)*X(0,0)*X(1,0)*X(2,1)
		 -M(0)*M(1)*V(0,2)*X(0,0)*X(1,1)*X(2,1)+M(1)*X(0,2)*X(1,1)*X(2,1)
		 -M(0)*M(2)*V(1,2)*X(0,0)*X(1,0)*X(2,2)+M(2)*X(0,2)*X(1,0)*X(2,2)
		 +M(1)*M(2)*V(0,2)*X(0,1)*X(1,1)*X(2,2)
		 -M(1)*M(2)*V(0,1)*X(0,2)*X(1,1)*X(2,2)
		 +M(0)*M(2)*V(0,1)*X(0,0)*X(1,2)*X(2,2)-M(2)*X(0,1)*X(1,2)*X(2,2)
		 -M(2)*V(0,2)*X(2,2)*Z(1,0)+M(2)*V(1,2)*X(2,2)*Z(1,0)
		 +M(1)*V(0,1)*X(1,1)*Z(2,0)-M(1)*V(1,2)*X(1,1)*Z(2,0)
		 -M(0)*V(0,1)*X(0,0)*Z(2,1)+M(0)*V(0,2)*X(0,0)*Z(2,1)
		 -X(1,0)*Z(2,0)+X(1,2)*Z(2,0)
		 +X(2,0)*Z(1,0)-X(2,1)*Z(1,0)
		 +X(0,1)*Z(2,1)-X(0,2)*Z(2,1));*/
}
