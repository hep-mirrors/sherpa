#ifndef MODEL_Interaction_Models_Interaction_Model_EW_Grav_H
#define MODEL_Interaction_Models_Interaction_Model_EW_Grav_H

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_EW_Grav : public Interaction_Model_Base {
    ATOOLS::Kabbala    g1,g2,sintW,costW,kap,om,vev;
    ATOOLS::Kabbala    PL,PR,M_I,root2,num2,num4,num15;
    ATOOLS::Kabbala K_CKM(short int i,short int j);  
    ATOOLS::Kabbala conj_K_CKM(short int i,short int j);
 public:
    Interaction_Model_EW_Grav(MODEL::Model_Base *,std::string,std::string);

    void c_FFT(std::vector<Single_Vertex>&,int&);
    void c_VVT(std::vector<Single_Vertex>&,int&);
    void c_SST(std::vector<Single_Vertex>&,int&);
    void c_VVVT(std::vector<Single_Vertex>&,int&);
    void c_FFVT(std::vector<Single_Vertex>&,int&);
    void c_SSST(std::vector<Single_Vertex>&,int&);

    void c_FFV(std::vector<Single_Vertex>&,int &)  {}
    void c_FFS(std::vector<Single_Vertex>&,int &)  {}
    void c_VVV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &)  {}
    void c_VVVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSSS(std::vector<Single_Vertex>&,int &) {}
  };
}

#endif





