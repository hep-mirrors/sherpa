// modul for parameters

#ifndef Model_MSSM_h
#define Model_MSSM_h

#include "Model_SM.H"
#include "Model_Higgs.H"
#include "Model_Chi.H"
#include "Model_sLeptons.H"
#include "Model_sLep_Chi.H"
#include "Model_SQCD.H"
#include "Model_SQCD_Chi.H"

namespace AMEGIC {  
  class Model_MSSM : public Model {
    Model_SM       mosm;
    Model_Higgs    moHiggs;
    Model_Chi      moch;
    Model_sLeptons mosL;
    Model_sLep_Chi mosLChi;
    Model_SQCD     mosq;
    Model_SQCD_Chi mosqChi;
  public:
    void Init();
    void c_FFV(Single_Vertex*,int&);
    void c_VVV(Single_Vertex*,int&);
    void c_FFS(Single_Vertex*,int&);
    void c_VVS(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);
    void c_SSV(Single_Vertex*,int&);
    
    inline double Aqcd(double);
    inline double Aqcd();
    inline double Aqed(double t);
    inline double Aqed();
    inline double SinTW();
    inline double CosTW();
    inline double TanB();
    inline double CosA();
    inline double SinA();
    //necessary for Amegic.C
    inline Spectrum_sUpquarks* Get_Spectrum_sUpquarks()     {return mosq.Get_Spectrum_sUpquarks();}
    inline Spectrum_sDownquarks* Get_Spectrum_sDownquarks() {return mosq.Get_Spectrum_sDownquarks();}
    inline Spectrum_sLeptons* Get_Spectrum_sLeptons()       {return mosL.Get_Spectrum_sLeptons();}
    inline Spectrum_Neutralinos* Get_Spectrum_Neutralinos() {return moch.Get_Spectrum_Neutralinos();}
    inline Spectrum_Charginos* Get_Spectrum_Charginos()     {return moch.Get_Spectrum_Charginos();}
    inline Spectrum_Higgs* Get_Spectrum_Higgs()             {return moHiggs.Get_Spectrum_Higgs();}
  };
}

#endif











