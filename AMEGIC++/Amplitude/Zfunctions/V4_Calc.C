#include "AMEGIC++/Amplitude/Zfunctions/Calculator.H"
#include "AMEGIC++/String/String_Generator.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

DEFINE_ZFCALC_GETTER(V4_Calc,"V4","v4 calculator")

V4_Calc::V4_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS) 
{ 
  type="V4";
  ncoupl=9;narg=8;pn=4;
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gamma",LF_Key()));
  lorentzlist.push_back(LF_Getter::GetObject("Gauge4",LF_Key()));
  for (short int i=0;i<4;i++) lorentzlist[i]->SetParticleArg(i);
  lorentzlist[4]->SetParticleArg(0,1,2,3);     
}

Kabbala V4_Calc::Massless()
{ return 2*Z(0,1)*Z(2,3)-Z(0,2)*Z(1,3)-Z(0,3)*Z(1,2);}


Kabbala V4_Calc::Do() 
{
  Kabbala factor = sgen->GetEnumber(coupl[8]);
 
  if (IsZero(M(0)) &&
      IsZero(M(1)) &&
      IsZero(M(2)) &&
      IsZero(M(3))) {
    return factor*Massless();
  }
  if (!IsZero(X(0,0)*M(0)) || 
      !IsZero(X(1,1)*M(1)) || 
      !IsZero(X(2,2)*M(2)) || 
      !IsZero(X(3,3)*M(3))) {
    std::cerr<<"Error in V4_Calc::Do(): not cutted massive vertex!"<<std::endl;
    abort();
  }
  return factor*Massless();

  /*  Kabbala M_1 = M(1)*X(1,1);
  Kabbala M_2 = M(2)*X(2,2);
  Kabbala M_3 = M(3)*X(3,3);
  Kabbala mt=
    M_3*(X(1,3)*Z(2,0) + X(0,3)*Z(2,1) - 2*Z(1,0)*X(2, 3)) 
    + M_2*(X(1,2)*Z(3,0) + X(0,2)*Z(3,1) 
	   - 2*Z(1,0)*(X(3,2) - M_3*V(2,3))
	   - M_3*(X(0,3)*X(1,2) + X(0,2)*X(1,3)))
    + M_1*(X(3,1)*Z(2,0) + X(2,1)*Z(3,0) - 2*X(0,1)*Z(3,2)
	   - M_3*(X(0,3)*X(2,1) - 2*X(0,1)*X(2,3) + V(1,3)*Z(2,0))
	   + M_2*(X(0,2)*(M_3*V(1,3) - X(3,1)) + 
		  2*X(0,1)*(X(3,2) - M_3*V(2,3)) + 
		  V(1,2)*(M_3*X(0,3) - Z(3,0)))) 
    + M(0)*X(0,0)*(X(3,0)*Z(2,1) + X(2,0)*Z(3,1) - 2*X(1,0)*Z(3,2)
		   - M_3*(X(1,3)*X(2,0) - 2*X(1,0)*X(2,3) + V(0,3)*Z(2,1))
		   + M_2*(X(1,2)*(M_3*V(0,3)- X(3,0)) +
			  2*X(1,0)*(X(3,2) - M_3*V(2,3)) +   
			  V(0,2)*(M_3*X(1,3) - Z(3,1)))    
		   + M_1*(2*V(0,1)*Z(3,2) - X(2,0)*X(3,1) 
			  + M_3*(V(1,3)*X(2,0) - 2*V(0,1)*X(2,3)) 
			  + X(2,1)*(M_3*V(0,3) - X(3,0))
			  + M_2*(V(1,2)*(X(3,0) - M_3*V(0,3)) + 
				 V(0,2)*(X(3,1) - M_3*V(1,3)) + 
				 2*V(0,1)*(M_3*V(2,3) - X(3,2)))))
				 - Z(2,0)*Z(3,1) - Z(2,1)*Z(3,0)+ 2*Z(1,0)*Z(3,2);*/


      /* 2*M(1)*M(2)*X(0,1)*X(1,1)*X(2,2)*X(3,2) 
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
					 + 2*V(0,1)*Z(3,2)));*/
  // return factor*mt;
}














