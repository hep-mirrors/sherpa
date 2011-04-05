#ifndef MODEL_Interaction_Models_Interaction_Model_SM_EHC_H
#define MODEL_Interaction_Models_Interaction_Model_SM_EHC_H

#include "MODEL/Interaction_Models/Interaction_Model_SM.H"

namespace MODEL {  
  class Interaction_Model_SM_EHC : public Interaction_Model_Base {
    Interaction_Model_SM        * p_mosm;
    ATOOLS::Kabbala ghpp,ghgg,g1,g2,g3,sintW,costW,PL,PR,M_I,root2,vev;
  public:
    Interaction_Model_SM_EHC(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_FFT(std::vector<Single_Vertex>&,int &)  {}
    void c_VVT(std::vector<Single_Vertex>&,int &)  {}
    void c_SST(std::vector<Single_Vertex>&,int &)  {}
    void c_VVVT(std::vector<Single_Vertex>&,int &)  {}
    void c_FFVT(std::vector<Single_Vertex>&,int &)  {}
    void c_SSST(std::vector<Single_Vertex>&,int &)  {}
    ~Interaction_Model_SM_EHC();
};
}
#endif
