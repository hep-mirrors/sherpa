#ifndef Interaction_Model_MSSM_LQQ_h
#define Interaction_Model_MSSM_LQQ_h

#include "Interaction_Model_MSSM.H"
#include "Kabbala.H"

namespace MODEL {  
  //check the summing Convention !!!!
  class Interaction_Model_MSSM_LQQ : public Interaction_Model_Base {
    ATOOLS::Kabbala PL,PR;
    ATOOLS::Kabbala M_I,root2,invroot2,K_zero,num_half,num_1,num_2,num_3,num_4,num_6;
    
    ATOOLS::Kabbala K_CKM(short int i,short int j);
    ATOOLS::Kabbala conj_K_CKM(short int i,short int j);
    ATOOLS::Kabbala K_LQQ(short int l,short int u,short int d);
    ATOOLS::Kabbala conj_K_LQQ(short int l,short int u,short int d);
    ATOOLS::Kabbala K_Z_D(short int i,short int j);
    ATOOLS::Kabbala K_Z_U(short int i,short int j);

    int gen_sUp(ATOOLS::Flavour fl);
    int gen_sDown(ATOOLS::Flavour fl);

    Interaction_Model_MSSM * p_momssm;
  public:
    Interaction_Model_MSSM_LQQ(MODEL::Model_Base *,std::string,std::string);
 
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);

    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_VVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_VVVT(std::vector<Single_Vertex>&,int&) {} 
    void c_FFVT(std::vector<Single_Vertex>&,int&) {} 
    void c_SSST(std::vector<Single_Vertex>&,int&) {} 
 };
}

#endif


