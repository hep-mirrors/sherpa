// modul for parameters

#ifndef Model_EE_QCD_h
#define Model_EE_QCD_h

#include "Model_QCD.H"
#include "Kabbala.H"
#include "Spectrum_EW.H"
#include "Couplings_EW.H"

namespace AMEGIC {  
  class Model_EE_QCD : public Model {
    Model_QCD     moqcd;
    Spectrum_EW   SpEW;
    Couplings_EW  Cpl;
    Kabbala    g1,g2,PL,PR,M_I;
  public:
    void Init();
    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);

    inline double Aqcd(double);    
    inline double Aqcd();
    inline double Aqed(double t);
    inline double Aqed();
    inline double SinTW()         {return Cpl.SinThetaW();}
    inline double CosTW()         {return Cpl.CosThetaW();}
    inline Kabbala K_sinTW()      {return Kabbala(std::string("sin\\theta_W"),SinTW());}
    inline Kabbala K_cosTW()      {return Kabbala(std::string("cos\\theta_W"),CosTW());}    
  };
}

#endif
