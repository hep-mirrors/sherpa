// modul for parameters

#ifndef Model_EW_h
#define Model_EW_h

#include "Model.H"
#include "Kabbala.H"
#include "Spectrum_EW.H"
#include "Couplings_EW.H"

namespace AMEGIC {  
  class Model_EW : public Model {
    Spectrum_EW   SpEW;
    Couplings_EW  CplEW;
    Kabbala    g1,g2,PL,PR,M_I,root2,vev;
  public:
    void Init();
    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_FFS(Single_Vertex*,int&);
    void c_VVS(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);
    void c_SSVV(Single_Vertex*,int&);
    void c_SSSS(Single_Vertex*,int&);

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





