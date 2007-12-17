#include "Calculator.H"
#include "String_Generator.H"
#include "Message.H"

using namespace AMEGIC;
using namespace ATOOLS;
using namespace MODEL;

AnomalousV3_Calc::AnomalousV3_Calc(Virtual_String_Generator* _sgen,Basic_Sfuncs* _BS) : 
  Basic_Func(_sgen,_BS), 
  Zfunc_Calc(_sgen,_BS),
  Basic_Zfunc(_sgen,_BS), 
  Basic_Xfunc(_sgen,_BS), 
  Basic_Mfunc(_sgen,_BS), 
  Basic_Vfunc(_sgen,_BS),
  Basic_Epsilonfunc(_sgen,_BS)
{ 
  type     = zl::AV3;
  ncoupl=14;narg=6;pn=3;
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::Gamma));
  lorentzlist.push_back(Lorentz_Function(lf::AGauge3));
  for (short int i=0;i<3;i++) lorentzlist[i].SetParticleArg(i);
  lorentzlist[3].SetParticleArg(0,1,2);     
}

Kabbala AnomalousV3_Calc::Do() 
{
  if (!IsZero(X(0,0)) || 
      !IsZero(X(1,1)) || 
      !IsZero(X(2,2))) {
    std::cerr<<"Error in AnomalousV3_Calc::Do(): not cutted vertex!"<<std::endl;
    abort();
  }

  Kabbala f1 = sgen->GetEnumber(coupl[6]);
  Kabbala f2 = sgen->GetEnumber(coupl[7]);
  Kabbala f3 = sgen->GetEnumber(coupl[8]);
  Kabbala f4 = sgen->GetEnumber(coupl[9]);
  Kabbala f5 = sgen->GetEnumber(coupl[10]);
  Kabbala f6 = sgen->GetEnumber(coupl[11]);
  Kabbala f7 = sgen->GetEnumber(coupl[12]);

  Kabbala t1 = Z(2,1)*(X(0,1)-X(0,2))+Z(2,0)*X(1,2)-Z(0,1)*X(2,1);
  Kabbala t2 = Z(0,1)*X(2,0)-Z(2,0)*X(1,0);
  Kabbala t3 = X(2,1)*X(0,2)*X(1,0)-X(0,1)*X(1,2)*X(2,0)
    +Z(0,1)*(X(2,0)*V(2,1)-X(2,1)*V(2,0))
    +Z(2,0)*(X(1,2)*V(1,0)-X(1,0)*V(1,2))
    +Z(1,2)*(X(0,1)*V(0,2)-X(0,2)*V(0,1));
  Kabbala t4 = Z(0,1)*X(2,0)+Z(2,0)*X(1,0);
  
  int sarg[6];
  sarg[0]=ps[0].numb;
  sarg[1]=ps[1].numb;
  sarg[2]=ps[2].numb;
  sarg[3]=BS->GetPolNumber(arg[0],arg[1],GetPMass(arg[0],arg[1]));
  sarg[4]=BS->GetPolNumber(arg[4],arg[5],GetPMass(arg[4],arg[5]));
  sarg[5]=BS->GetPolNumber(arg[8],arg[9],GetPMass(arg[8],arg[9]));

//   std::cout<<sarg[0]<<" "<<sarg[1]<<" "<<sarg[2]<<" "<<sarg[3]<<" "<<sarg[4]<<" "<<sarg[5]<<std::endl;

  Kabbala t5 = Epsilon(sarg[5],sarg[4],sarg[1],sarg[3],1)-
               Epsilon(sarg[5],sarg[4],sarg[2],sarg[3],1);
  Kabbala t6 = Epsilon(sarg[5],sarg[4],sarg[0],sarg[3],1);
  Kabbala t7 = 
     Epsilon(sarg[4],sarg[2],sarg[0],sarg[3],1)*X(2,1)+
     Epsilon(sarg[1],sarg[5],sarg[0],sarg[3],1)*X(1,2)-
     Epsilon(sarg[4],sarg[5],sarg[0],sarg[3],1)*V(1,2)-
     Epsilon(sarg[1],sarg[2],sarg[0],sarg[3],1)*Z(1,2);
//    std::cout<<f5.Value()<<" "<<f6.Value()<<" "<<f7.Value()<<std::endl;

  return f1*t1+f2*t2-f3*t3+f4*t4+t5*f5+t6*f6-t7*f7;

}
