// modul for parameters

#ifndef Model_SQCD_h
#define Model_SQCD_h

#include "Model_Higgs.H"
#include "Spectrum_sDownquarks.H"
#include "Spectrum_sUpquarks.H"
 
namespace AMEGIC {  
  class Model_SQCD : public virtual Model {
    Spectrum_sDownquarks   SpsD;
    Spectrum_sUpquarks     SpsU;
    Kabbala    g1,g2,g3,PL,PR,M_I,root2,inv_root2,K_zero,num_2,num_3 ;    
    
    Model_Higgs* moHiggs;
    int          ini_higgs;
  public:
    Model_SQCD() : Model() { moHiggs=NULL;ini_higgs=0; }
    ~Model_SQCD();
    void Init();    
    void Init(Model_Higgs*,Isajet*);
    void c_SSV(Single_Vertex*,int&);
    void c_SSS(Single_Vertex*,int&);
    void c_FFS(Single_Vertex*,int&);
    void c_FFV(Single_Vertex*,int&);
    void c_SSVV(Single_Vertex*,int&);

    inline Kabbala K_sinTW()                       {return moHiggs->K_sinTW();}
    inline Kabbala K_cosTW()                       {return moHiggs->K_cosTW();}  
    inline Kabbala K_v1()                          {return moHiggs->K_v1();}
    inline Kabbala K_v2()                          {return moHiggs->K_v2();}    
    inline Kabbala K_Z_H(short int i,short int j)  {return moHiggs->K_Z_H(i,j);}     
    inline Kabbala K_Z_R(short int i,short int j)  {return moHiggs->K_Z_R(i,j);}     
    inline Kabbala K_B_R(short int i)              {return moHiggs->K_B_R(i);}
    inline Kabbala K_CKM(short int i,short int j)  {return moHiggs->K_CKM(i,j);}
    inline Kabbala conj_K_CKM(short int i,short int j)  {return moHiggs->conj_K_CKM(i,j);}
    inline Kabbala K_h()                           {return Kabbala(std::string("h"),SpsU.h());}
    inline Kabbala conj_K_h()                     {return Kabbala(std::string("h^\\ti"),SpsU.h());}
    //necessary for Amegic.C
    inline Spectrum_sUpquarks* Get_Spectrum_sUpquarks() {return &SpsU;}
    inline Spectrum_sDownquarks* Get_Spectrum_sDownquarks() {return &SpsD;}
    Kabbala K_w_S(short int i,short int j);
    Kabbala K_u_S(short int i,short int j);
    Kabbala K_e_S(short int i,short int j);
    Kabbala K_d_S(short int i,short int j);
    Kabbala K_Z_U(short int i,short int j);
    Kabbala K_Z_D(short int i,short int j);
    int gen_sUp(APHYTOOLS::Flavour fl);
    int gen_sDown(APHYTOOLS::Flavour fl);
 };
}

#endif

