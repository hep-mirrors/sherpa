// modul for parameters

#ifndef Model_SQCD_Chi_h
#define Model_SQCD_Chi_h

#include "Model_SM.H"
#include "Model_Chi.H"
#include "Model_SQCD.H"

namespace AMEGIC {  
  class Model_SQCD_Chi : public virtual Model {
    Kabbala    g1,g2,PL,PR,M_I,root2,K_zero,num_2,num_3;    
    
    Model_SM*       mosm;
    Model_Chi*      moch;
    Model_SQCD*     mosq;
    int             ini_sm,ini_ch,ini_sq;
  public:
    Model_SQCD_Chi() : Model() { mosm=NULL;moch=NULL;mosq=NULL;ini_sm=ini_ch=ini_sq=0; }
    ~Model_SQCD_Chi();
    void Init();    
    void Init(Model_SM*,Model_Chi*,Model_SQCD*);
    void c_FFS(Single_Vertex*,int&);
    
    inline Kabbala K_sinTW()                       {return moch->K_sinTW();}
    inline Kabbala K_cosTW()                       {return moch->K_cosTW();}     
    inline Kabbala K_yuk_sign(APHYTOOLS::Flavour fl)          {return moch->K_yuk_sign(fl);}        
    inline Kabbala K_Z_PL(short int i,short int j) {return moch->K_Z_PL(i,j);}
    inline Kabbala K_Z_N(short int i,short int j)  {return moch->K_Z_N(i,j);}
    inline Kabbala K_Z_MI(short int i,short int j) {return moch->K_Z_MI(i,j);}    
    inline Kabbala K_Z_U(short int i,short int j)  {return mosq->K_Z_U(i,j);}
    inline Kabbala K_Z_D(short int i,short int j)  {return mosq->K_Z_D(i,j);}
    inline Kabbala K_CKM(short int i,short int j)  {return mosm->K_CKM(i,j);} 
    inline int gen_sUp(APHYTOOLS::Flavour fl)                 {return mosq->gen_sUp(fl);}
    inline int gen_sDown(APHYTOOLS::Flavour fl)               {return mosq->gen_sDown(fl);}
    inline Spectrum_sUpquarks* Get_Spectrum_sUpquarks()     {return mosq->Get_Spectrum_sUpquarks();}
    inline Spectrum_sDownquarks* Get_Spectrum_sDownquarks() {return mosq->Get_Spectrum_sDownquarks();}
    inline Spectrum_Neutralinos* Get_Spectrum_Neutralinos() {return moch->Get_Spectrum_Neutralinos();}
    inline Spectrum_Charginos* Get_Spectrum_Chargions()     {return moch->Get_Spectrum_Charginos();}    
    inline Kabbala K_CKM_conj(short int i,short int j);    
  };
}

#endif

