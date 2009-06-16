#ifndef MODEL_Interaction_Models_Interaction_Model_sLepton_sQuark_H
#define MODEL_Interaction_Models_Interaction_Model_sLepton_sQuark_H

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  //check the summing Convention !!!!
  class Interaction_Model_sLepton_sQuark : public Interaction_Model_Base {
    ATOOLS::Kabbala g1,g2,vev,v1,v2,sintW,costW,PL,PR,mu,conj_mu;
    ATOOLS::Kabbala M_I,root2,invroot2,K_zero,num_2,num_3,num_4,num_6;
    
    ATOOLS::Kabbala K_Z_Nu(short int i,short int j);
    ATOOLS::Kabbala K_Z_L(short int i,short int j);
    ATOOLS::Kabbala K_Z_U(short int i,short int j);
    ATOOLS::Kabbala K_Z_D(short int i,short int j);
    ATOOLS::Kabbala K_CKM(short int i,short int j);
    ATOOLS::Kabbala conj_K_CKM(short int i,short int j);
    ATOOLS::Kabbala K_l(short int i);
    ATOOLS::Kabbala K_u(short int i);
    ATOOLS::Kabbala K_d(short int i);
    ATOOLS::Kabbala K_yuk(ATOOLS::Flavour fl);           
    ATOOLS::Kabbala K_yuk_sign(ATOOLS::Flavour fl);          
    int gen_sLep(ATOOLS::Flavour fl);
    int gen_sUp(ATOOLS::Flavour fl);
    int gen_sDown(ATOOLS::Flavour fl);
  public:
    Interaction_Model_sLepton_sQuark(MODEL::Model_Base *,std::string,std::string);
 
    void c_FFV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVVV(std::vector<Single_Vertex>&,int &) {}    
    void c_FFS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSVV(std::vector<Single_Vertex>&,int &) {}
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


