// modul for parameters

#ifndef Model_EW_Grav_h
#define Model_EW_Grav_h

#include "Model.H"
#include "Kabbala.H"
#include "Spectrum_EW.H"
#include "Couplings_EW.H"
#include "Couplings_LED.H"

namespace AMEGIC {  
  class Model_EW_Grav : public Model {
    Spectrum_EW   SpEW;
    Couplings_EW  CplEW;
    Couplings_LED  CplLED;
    Kabbala    g1,g2,PL,PR,M_I,root2,vev;
    Kabbala    kap,om,num2,num4,num15;
  public:
    void Init();
    
    void c_FFT(Single_Vertex*,int&);
    void c_VVT(Single_Vertex*,int&);
    void c_SST(Single_Vertex*,int&);
    void c_VVVT(Single_Vertex*,int&);
    void c_FFVT(Single_Vertex*,int&);
    void c_SSST(Single_Vertex*,int&);
      
      inline double Aqed(double t);
      inline double Aqed();
      inline double SinTW()           {return CplEW.SinThetaW();}
      inline double CosTW()           {return CplEW.CosThetaW();}
      inline double v_0()             {return std::abs(vev.Value());}
      inline Spectrum_EW* Get_Spectrum_EW() {return &SpEW;}
      inline Kabbala K_sinTW()        {return Kabbala(std::string("sin\\theta_W"),SinTW());}
      inline Kabbala K_cosTW()        {return Kabbala(std::string("cos\\theta_W"),CosTW());}
      Complex CKM(short int i,short int j);
      Kabbala K_CKM(short int i,short int j);  
      Kabbala conj_K_CKM(short int i,short int j);
 };
}

#endif





