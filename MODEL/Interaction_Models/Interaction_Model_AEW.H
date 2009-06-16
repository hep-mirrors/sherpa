#ifndef MODEL_Interaction_Models_Interaction_Model_AEW_H
#define MODEL_Interaction_Models_Interaction_Model_AEW_H

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_AEW : public Interaction_Model_Base {
    ATOOLS::Kabbala g1,g2,sintW,costW,PL,PR,M_I,root2,vev,a4,a5;
    ATOOLS::Kabbala g1_p,lambda_p,kappa_p,g4_p,g5_p,lambdat_p,kappat_p;
    ATOOLS::Kabbala g1_Z,lambda_Z,kappa_Z,g4_Z,g5_Z,lambdat_Z,kappat_Z;
    ATOOLS::Kabbala f4_Z,f5_Z,f4_p,f5_p;
    ATOOLS::Kabbala h1_Z,h2_Z,h3_Z,h4_Z,h1_p,h2_p,h3_p,h4_p;
    ATOOLS::Kabbala K_CKM(short int i,short int j);  
    ATOOLS::Kabbala conj_K_CKM(short int i,short int j);
  public:
    Interaction_Model_AEW(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &) {}
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &) {}
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &) {}
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSSS(std::vector<Single_Vertex>&,int &) {}

    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_VVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_VVVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_FFVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SSST(std::vector<Single_Vertex>&,int&)  {} 
  };
}

#endif


