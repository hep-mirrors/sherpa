// modul for parameters

#ifndef Model_Higgs_h
#define Model_Higgs_h

#include "Model_SM.H"
#include "Kabbala.H"
#include "Spectrum_Higgs.H"

namespace AMEGIC {  
  class Model_Higgs : public Model {
    Spectrum_Higgs       SpHiggs;
    Model_SM*            mosm;
    Kabbala              g1,g2,PL,PR,M_I,root2,K_zero,num_2,num_3,num_4;    
    int ini_sm;
 public:
    Model_Higgs() : Model() { mosm=NULL;ini_sm=0; }
    ~Model_Higgs();    
    void Init();
    void Init(Model_SM*,Isajet*);
    void c_FFS(Single_Vertex*,int&);
    void c_VVS(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);
    void c_SSV(Single_Vertex*,int&);
    
    inline Kabbala K_sinTW()                     {return mosm->K_sinTW();}
    inline Kabbala K_cosTW()                     {return mosm->K_cosTW();}    
    inline double Aqed()                         {return mosm->Aqed();}
    inline double Aqcd()                         {return mosm->Aqcd();}
    inline double SinTW()                        {return mosm->SinTW();}
    inline double CosTW()                        {return mosm->CosTW();}    
    inline double TanB()                         {return SpHiggs.Z_H(0,0)/SpHiggs.Z_H(1,0);}
    inline double CosA()                         {return SpHiggs.Z_R(0,1);}
    inline double SinA()                         {return SpHiggs.Z_R(1,1);}    
    inline Kabbala K_CKM(short int i,short int j){return mosm->K_CKM(i,j);}           
    inline Kabbala conj_K_CKM(short int i,short int j){return mosm->conj_K_CKM(i,j);}           
    Kabbala K_yuk(APHYTOOLS::Flavour fl);           
    Kabbala K_yuk_sign(APHYTOOLS::Flavour fl);          
    Kabbala K_v0();     
    Kabbala K_v1();    
    Kabbala K_v2();    
    Kabbala K_Z_H(short int i,short int j);    
    Kabbala K_Z_R(short int i,short int j);    
    Kabbala K_A_M(short int i,short int j);    
    Kabbala K_A_H(short int i,short int j);        
    Kabbala K_A_P(short int i,short int j);    
    Kabbala K_B_R(short int i);

    inline Spectrum_Higgs* Get_Spectrum_Higgs() {return &SpHiggs;}   
  };
}

#endif

