#include "Calculator.H"
#include "String_Generator.H"

using namespace AMEGIC;

V4_Calc::V4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type     = zl::V4;
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gauge4));
  for (short int i=0;i<4;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[4].SetParticleArg(0,1,2,3);     
}

Kabbala V4_Calc::Massless()
{ return 2*Z(0,1)*Z(2,3)-Z(0,2)*Z(1,3)-Z(0,3)*Z(1,2);}


Kabbala V4_Calc::Do() 
{
  Kabbala factor = sgen->Get_Enumber(coupl[8]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3))) {
    return factor*Massless();
  }
  

  return factor*(
     2*M(1)*M(2)*X(0,1)*X(1,1)*X(2,2)*X(3,2) 
    -M(1)*M(3)*X(0,3)*X(1,1)*X(2,1)*X(3,3) 
    -2*M(1)*M(2)*M(3)*V(2,3)*X(0,1)*X(1,1)*X(2,2)*X(3,3) 
    +M(1)*M(2)*M(3)*V(1,2)*X(0,3)*X(1,1)*X(2,2)*X(3,3) 
    -M(2)*M(3)*X(0,3)*X(1,2)*X(2,2)*X(3,3) 
    +2*M(1)*M(3)*X(0,1)*X(1,1)*X(2,3)*X(3,3) 
    -2*M(2)*X(2,2)*X(3,2)*Z(1,0) 
    +2*M(2)*M(3)*V(2,3)*X(2,2)*X(3,3)*Z(1,0) 
    -2*M(3)*X(2,3)*X(3,3)*Z(1,0) 
    +M(1)*X(1,1)*X(3,1)*Z(2,0) 
    -M(1)*M(3)*V(1,3)*X(1,1)*X(3,3)*Z(2,0) 
    +M(3)*X(1,3)*X(3,3)*Z(2,0) 
    +M(3)*X(0,3)*X(3,3)*Z(2,1) 
    +M(1)*X(1,1)*X(2,1)*Z(3,0) 
    -M(1)*M(2)*V(1,2)*X(1,1)*X(2,2)*Z(3,0) 
    +M(2)*X(1,2)*X(2,2)*Z(3,0) 
    -Z(2,1)*Z(3,0) 
    -Z(2,0)*Z(3,1) 
    +M(2)*X(0,2)*X(2,2)*(-(M(3)*X(1,3)*X(3,3)) 
			 +X(1,1)*(-(M(1)*X(3,1)) 
				  + M(1)*M(3)*V(1,3)*X(3,3)) 
			 + Z(3,1)) 
    -2*M(1)*X(0,1)*X(1,1)*Z(3,2) 
    +2*Z(1,0)*Z(3,2) 
    +M(0)*X(0,0)*(2*M(2)*X(1,0)*X(2,2)*X(3,2) 
		  -M(3)*X(1,3)*X(2,0)*X(3,3) 
		  -2*M(2)*M(3)*V(2,3)*X(1,0)*X(2,2)*X(3,3) 
		  +M(2)*M(3)*V(0,2)*X(1,3)*X(2,2)*X(3,3) 
		  +2*M(3)*X(1,0)*X(2,3)*X(3,3) 
		  +M(2)*X(1,2)*X(2,2)*(-X(3,0) 
				       + M(3)*V(0,3)*X(3,3)) 
		  +X(3,0)*Z(2,1) 
		  -M(3)*V(0,3)*X(3,3)*Z(2,1) 
		  +X(2,0)*Z(3,1) 
		  -M(2)*V(0,2)*X(2,2)*Z(3,1) 
		  -2*X(1,0)*Z(3,2) 
		  +M(1)*X(1,1)*(-(X(2,0)*X(3,1)) 
				+M(2)*V(0,2)*X(2,2)*X(3,1) 
				-2*M(2)*V(0,1)*X(2,2)*X(3,2) 
				+M(3)*V(1,3)*X(2,0)*X(3,3) 
				-M(2)*M(3)*V(0,2)*V(1,3)*X(2,2)*X(3,3) 
				+2*M(2)*M(3)*V(0,1)*V(2,3)*X(2,2)*X(3,3) 
				-2*M(3)*V(0,1)*X(2,3)*X(3,3) 
				+M(2)*V(1,2)*X(2,2)*(X(3,0) 
						     - M(3)*V(0,3)*X(3,3)) 
				+X(2,1)*(-X(3,0) 
					 + M(3)*V(0,3)*X(3,3)) 
				+ 2*V(0,1)*Z(3,2))));
}














