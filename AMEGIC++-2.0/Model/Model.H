// modul for parameters

#ifndef Model_h
#define Model_h

#include "Vertex.H"
#include "Isajet/Isajet.H"
#include "Spectrum_EW.H"
#include "Spectrum_Higgs.H"
#include "Spectrum_sUpquarks.H"
#include "Spectrum_sDownquarks.H"
#include "Spectrum_sLeptons.H"
#include "Spectrum_Neutralinos.H"
#include "Spectrum_Charginos.H"

namespace AMEGIC {

  struct Single_Vertex; 
  class Model {
    Vertex* v;
  protected:  
    Isajet* isa;
    int     ini_isa;
  public:
    Model() { v=NULL;isa=NULL;ini_isa = 0;} 
    Vertex* Get_Vertex() {return v;}
    void    Init_Vertex();
    void    Init_Isajet() {
      if (isa!=NULL) delete isa;
      ini_isa = 1;
      isa = new Isajet;
    } 

    virtual ~Model() {
      if (v   != 0) { delete v; v = 0; }
      if (isa != 0) { delete isa; isa = 0; }
    }
    virtual void Init();
    virtual void  c_FFV(Single_Vertex*,int&)  { }
    virtual void  c_FFS(Single_Vertex*,int&)  { }
    virtual void  c_VVV(Single_Vertex*,int&)  { }
    virtual void  c_SSV(Single_Vertex*,int&)  { }
    virtual void  c_VVS(Single_Vertex*,int&)  { }
    virtual void  c_SSS(Single_Vertex*,int&)  { }
    virtual void  c_VVVV(Single_Vertex*,int&) { }
    virtual void  c_SSVV(Single_Vertex*,int&) { }
    virtual void  c_SSSS(Single_Vertex*,int&) { }
    virtual double SinTW();
    virtual double CosTW();
    virtual double Aqed();
    virtual double Aqed(double);
    virtual double Aqcd();
    virtual double Aqcd(double);
    virtual double TanB();
    virtual double CosA();
    virtual double SinA();
    virtual Spectrum_EW*  Get_Spectrum_EW();
    virtual Spectrum_Higgs*  Get_Spectrum_Higgs();
    virtual Spectrum_sUpquarks*  Get_Spectrum_sUpquarks();
    virtual Spectrum_sDownquarks*  Get_Spectrum_sDownquarks();
    virtual Spectrum_sLeptons*  Get_Spectrum_sLeptons();
    virtual Spectrum_Neutralinos*  Get_Spectrum_Neutralinos();
    virtual Spectrum_Charginos*  Get_Spectrum_Charginos();
  };
}

//#include "Vertex.H"

extern AMEGIC::Model* mo;

#endif







