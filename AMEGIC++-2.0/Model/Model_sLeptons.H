// modul for parameters

#ifndef Model_sLeptons_h
#define Model_sLeptons_h

#include "Model_Higgs.H"
#include "Spectrum_sNeutrinos.H"
#include "Spectrum_sLeptons.H"

namespace AMEGIC {  
  class Model_sLeptons : public virtual Model {
    Spectrum_sNeutrinos  SpsN;
    Spectrum_sLeptons    SpsL;
    Kabbala              g1,g2,PL,PR,M_I,root2,K_zero,num_2,num_4;    
    
    Model_Higgs* moHiggs;
    int          ini_higgs;
  public:
    Model_sLeptons() : Model() { moHiggs=NULL;ini_higgs=0; }
    ~Model_sLeptons();
    void Init();    
    void Init(Model_Higgs*,Isajet*);
    void c_SSS(Single_Vertex*,int&);
    void c_SSV(Single_Vertex*,int&);

    inline Kabbala K_sinTW()                        {return moHiggs->K_sinTW();}
    inline Kabbala K_cosTW()                        {return moHiggs->K_cosTW();}    
    inline Kabbala K_v1()                           {return moHiggs->K_v1();}    
    inline Kabbala K_v2()                           {return moHiggs->K_v2();}            
    inline Kabbala K_Z_R(short int i,short int j)   {return moHiggs->K_Z_R(i,j);}   
    inline Kabbala K_Z_H(short int i,short int j)   {return moHiggs->K_Z_H(i,j);} 
    inline Kabbala K_yuk(APHYTOOLS::Flavour fl)                {return moHiggs->K_yuk(fl);} 
    inline Kabbala K_yuk_sign(APHYTOOLS::Flavour fl)           {return moHiggs->K_yuk_sign(fl);}            
    inline Kabbala K_B_R(short int i)               {return moHiggs->K_B_R(i);}
    inline Kabbala K_h()                            {return Kabbala(std::string("h"),SpsL.h());}
    //necessary for Amegic.C
    inline Spectrum_sLeptons* Get_Spectrum_sLeptons(){return &SpsL;}
    inline Kabbala K_l_S(short int i,short int j);
    inline Kabbala K_k_S(short int i,short int j);
    Kabbala K_Z_nue(short int i,short int j);    
    Kabbala K_Z_L(short int i,short int j);    
    int gen_sLep(APHYTOOLS::Flavour fl);
  };
}

#endif


