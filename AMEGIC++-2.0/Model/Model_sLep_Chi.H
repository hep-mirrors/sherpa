// modul for parameters

#ifndef Model_sLep_Chi_h
#define Model_sLep_Chi_h

#include "Model_SM.H"
#include "Model_Chi.H"
#include "Model_sLeptons.H"

namespace AMEGIC {  
  class Model_sLep_Chi : public virtual Model {
    Kabbala    g1,g2,PL,PR,M_I,root2,K_zero;    
    
    Model_SM*       mosm;
    Model_Chi*      moch;
    Model_sLeptons* mosL;
    int             ini_sm,ini_ch,ini_sL;
  public:
    Model_sLep_Chi() : Model() { mosm=NULL;moch=NULL;mosL=NULL;ini_sm=ini_ch=ini_sL=0; }
    ~Model_sLep_Chi();
    void Init();    
    void Init(Model_SM*,Model_Chi*,Model_sLeptons*);
    void c_FFS(Single_Vertex*,int&);
    
    inline Kabbala K_sinTW()                                {return moch->K_sinTW();}
    inline Kabbala K_cosTW()                                {return moch->K_cosTW();}    
    inline Kabbala K_Z_PL(short int i,short int j)          {return moch->K_Z_PL(i,j);}    
    inline Kabbala K_Z_MI(short int i,short int j)          {return moch->K_Z_MI(i,j);}    
    inline Kabbala K_Z_N(short int i,short int j)           {return moch->K_Z_N(i,j);}   
    inline Kabbala K_Z_N_com(short int i,short int j)       {return moch->K_Z_N_com(i,j);}   
    inline Kabbala K_Z_N_com_conj(short int i,short int j)  {return moch->K_Z_N_com_conj(i,j);}   
    inline Kabbala K_Z_nue(short int i,short int j)         {return mosL->K_Z_nue(i,j);}    
    inline Kabbala K_Z_L(short int i,short int j)           {return mosL->K_Z_L(i,j);}    
    inline Kabbala K_yuk_sign(APHYTOOLS::Flavour fl)                   {return mosL->K_yuk_sign(fl);}
    inline int gen_sLep(APHYTOOLS::Flavour fl)                         {return mosL->gen_sLep(fl);} 
    inline Spectrum_sLeptons* Get_Spectrum_sLeptons()       {return mosL->Get_Spectrum_sLeptons();}
    inline Spectrum_Neutralinos* Get_Spectrum_Neutralinos() {return moch->Get_Spectrum_Neutralinos();}
    inline Spectrum_Charginos* Get_Spectrum_Charginos()     {return moch->Get_Spectrum_Charginos();}
  }; 
}
 
#endif


