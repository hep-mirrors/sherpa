#ifndef Interaction_Model_Inos_h
#define Interaction_Model_Inos_h

#include "Interaction_Model_Base.H"
#include "Kabbala.H"


namespace AMEGIC {  
  class Interaction_Model_Inos : public Interaction_Model_Base {
    Kabbala g1,g2,sintW,costW,PL,PR,M_I,root2,K_zero,num_2,num_4;
    
    Kabbala K_CKM(short int i,short int j);
    Kabbala conj_K_CKM(short int i,short int j);
    Kabbala K_Z_H(short int i,short int j);
    Kabbala K_Z_R(short int i,short int j);
    Kabbala K_yuk(APHYTOOLS::Flavour fl);
    Kabbala K_yuk_sign(APHYTOOLS::Flavour fl);
    Kabbala K_Z_PL(short int i,short int j);    
    Kabbala K_Z_MI(short int i,short int j);    
    Kabbala K_Z_N(short int i,short int j);    
    Kabbala K_Z_N_com(short int i,short int j);       
    Kabbala K_Z_N_com_conj(short int i,short int j);       
  public:
    Interaction_Model_Inos(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(Single_Vertex *,int &);
    void c_VVV(Single_Vertex *,int &)  {}
    void c_VVVV(Single_Vertex *,int &) {}
    void c_FFS(Single_Vertex *,int &);
    void c_SSV(Single_Vertex *,int &)  {}
    void c_VVS(Single_Vertex *,int &)  {}
    void c_SSS(Single_Vertex *,int &)  {}
    void c_SSVV(Single_Vertex *,int &) {}
    void c_SSSS(Single_Vertex *,int &) {}
    
    void c_FFT(Single_Vertex  *,int&)  {}
    void c_VVT(Single_Vertex  *,int&)  {} 
    void c_SST(Single_Vertex  *,int&)  {} 
    void c_VVVT(Single_Vertex *,int&)  {} 
    void c_FFVT(Single_Vertex *,int&)  {} 
    void c_SSST(Single_Vertex *,int&)  {} 
  };
}

#endif

