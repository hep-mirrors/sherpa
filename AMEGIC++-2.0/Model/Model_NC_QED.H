// modul for parameters

#ifndef Model_NC_QED_h
#define Model_NC_QED_h

#include "Model.H"
#include "Kabbala.H"
#include "Spectrum_EW.H"
#include "Couplings_EW.H"

namespace AMEGIC {  
  class Model_NC_QED : public Model {
    Spectrum_EW   SpEW;
    Couplings_EW  Cpl;
    Kabbala       g1,PL,PR,M_I,num2;
  public:
    void Init();
    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);

    inline double Aqed(double t);
    inline double Aqed();
  };
}

#endif
