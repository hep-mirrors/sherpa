// modul for parameters

#ifndef Model_SM_h
#define Model_SM_h

#include "Model_QCD.H"
#include "Model_EW.H"

namespace AMEGIC {  
  class Model_SM : public Model {
    Model_QCD moqcd;
    Model_EW  moew;
  public:
    void Init();
    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);
    void c_FFS(Single_Vertex*,int&);
    void c_VVS(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);

    inline double Aqcd(double );    
    inline double Aqcd();           
    inline double Aqed(double t);  
    inline double Aqed();          
    inline double SinTW()                             {return moew.SinTW();}
    inline double CosTW()                             {return moew.CosTW();}
    //    inline Spectrum_QCD* Get_Spectrum_QCD()           {return moqcd.Get_Spectrum_QCD();}
    inline Spectrum_EW*  Get_Spectrum_EW()            {return moew.Get_Spectrum_EW();}
    inline Kabbala K_sinTW()                          {return moew.K_sinTW();}
    inline Kabbala K_cosTW()                          {return moew.K_cosTW();}    
    inline Complex CKM(short int i,short int j)       {return moew.CKM(i,j);}
    inline Kabbala K_CKM(short int i,short int j)     {return moew.K_CKM(i,j);}
    inline Kabbala conj_K_CKM(short int i,short int j){return moew.K_CKM(i,j);}
 };
}

#endif
