// modul for parameters

#ifndef Model_THDM_h
#define Model_THDM_h

#include "Model_SM.H"
#include "Model_Higgs.H"

namespace AMEGIC {  
  class Model_THDM : public Model {
    Model_SM    mosm;
    Model_Higgs moHiggs;
  public:
    void Init();
    
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_FFS(Single_Vertex*,int&);
    void c_VVS(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);
    void c_SSV(Single_Vertex*,int&);
    void c_SSVV(Single_Vertex*,int&);
    void c_VVVV(Single_Vertex*,int&);
    
    inline double Aqcd(double);
    inline double Aqcd();
    inline double Aqed(double t);
    inline double Aqed();
    inline double SinTW()          {return mosm.SinTW();}
    inline double CosTW()          {return mosm.CosTW();}
    inline double TanB()           {return moHiggs.TanB();}
    inline double CosA()           {return moHiggs.CosA();}
    inline double SinA()           {return moHiggs.SinA();}
  };
}

#endif











