// modul for parameters

#ifndef Model_QCD_h
#define Model_QCD_h

#include "Model.H"
#include "Kabbala.H"
//#include "Spectrum_QCD.H"
// #include "Couplings_QCD.H"

namespace AMEGIC {  
  class Model_QCD : public virtual Model {
    // Spectrum_QCD   SpQCD;
    // Couplings_QCD  Cpl;
    Kabbala        g3,PL,PR,M_I;    
  public:
    void Init();    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);
    inline double Aqcd(double);
    inline double Aqcd();
    //    inline Spectrum_QCD* Get_Spectrum_QCD() {return &SpQCD;}
  };
}

#endif


