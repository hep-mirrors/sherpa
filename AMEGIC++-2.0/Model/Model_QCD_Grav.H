// modul for parameters

#ifndef Model_QCD_Grav_h
#define Model_QCD_Grav_h

#include "Model.H"
#include "Kabbala.H"
#include "Couplings_LED.H"
// #include "Couplings_QCD.H"

namespace AMEGIC {  
  class Model_QCD_Grav : public virtual Model {
    Couplings_LED  CplLED;
    // Couplings_QCD  Cpl;
    Kabbala        g3,PL,PR,M_I;    
    Kabbala    kap,om,num2,num4,num15;
  public:
    void Init();    
    void c_FFVT(Single_Vertex*,int&);
    void c_VVT(Single_Vertex*,int&);
    void c_VVVT(Single_Vertex*,int&);
    inline double Aqcd(double);
    inline double Aqcd();
  };
}

#endif


