#ifndef Interaction_Model_Higgs_h
#define Interaction_Model_Higgs_h

#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace AMEGIC {  
  class Interaction_Model_Higgs : public Interaction_Model_Base {
    Kabbala g1,g2,sintW,costW,PL,PR,M_I,root2,vev,v1,v2,K_zero,num_2,num_3,num_4;    

    Kabbala K_CKM(short int i,short int j);  
    Kabbala conj_K_CKM(short int i,short int j);
    Kabbala K_yuk(ATOOLS::Flavour fl);           
    Kabbala K_yuk_sign(ATOOLS::Flavour fl);          
    Kabbala K_Z_H(short int i,short int j);    
    Kabbala K_Z_R(short int i,short int j);    
    Kabbala K_A_M(short int i,short int j);    
    Kabbala K_A_H(short int i,short int j);        
    Kabbala K_A_P(short int i,short int j);    
    Kabbala K_A_R(short int i,short int j);    
    Kabbala K_B_R(short int i);
  public:
    Interaction_Model_Higgs(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &)  {}
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &);
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);
  
    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_VVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_VVVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_FFVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SSST(std::vector<Single_Vertex>&,int&)  {} 
  };
}

#endif


