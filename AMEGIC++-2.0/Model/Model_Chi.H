// modul for parameters

#ifndef Model_Chi_h
#define Model_Chi_h

#include "Model_Higgs.H"
#include "Spectrum_Charginos.H"
#include "Spectrum_Neutralinos.H"
 
namespace AMEGIC {  
  class Model_Chi : public virtual Model {
    Spectrum_Charginos   SpCh;
    Spectrum_Neutralinos SpNe;
    Kabbala    g1,g2,PL,PR,M_I,root2,K_zero,num_2,num_4;    
    
    Model_Higgs* moHiggs;
    int          ini_higgs;
  public:
    Model_Chi() : Model() { moHiggs=NULL;ini_higgs=0; }
    ~Model_Chi();
    void Init();    
    void Init(Model_Higgs*,Isajet*);
    void c_FFS(Single_Vertex*,int&);
    void c_FFV(Single_Vertex*,int&);
    
    inline Kabbala K_sinTW()                       {return moHiggs->K_sinTW();}
    inline Kabbala K_cosTW()                       {return moHiggs->K_cosTW();}    
    inline Kabbala K_Z_H(short int i,short int j)  {return moHiggs->K_Z_H(i,j);}     
    inline Kabbala K_Z_R(short int i,short int j)  {return moHiggs->K_Z_R(i,j);}     
    inline Kabbala K_yuk_sign(APHYTOOLS::Flavour fl)          {return moHiggs->K_yuk_sign(fl);}    
    inline Spectrum_Neutralinos*  Get_Spectrum_Neutralinos()  {return &SpNe;}
    inline Spectrum_Charginos*  Get_Spectrum_Charginos()      {return &SpCh;}
    Kabbala K_Z_PL(short int i,short int j);    
    Kabbala K_Z_MI(short int i,short int j);    
    Kabbala K_Z_N(short int i,short int j);    
    Kabbala K_Z_N_com(short int i,short int j);       
    Kabbala K_Z_N_com_conj(short int i,short int j);       
 };
}

#endif

