#ifndef MODEL_Interaction_Models_Interaction_Model_QCD_Grav_H
#define MODEL_Interaction_Models_Interaction_Model_QCD_Grav_H

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_QCD_Grav : public Interaction_Model_Base {
    ATOOLS::Kabbala        g3,PL,PR,M_I;    
    ATOOLS::Kabbala        kap,om,num2,num4,num15;
  public:
    Interaction_Model_QCD_Grav(MODEL::Model_Base *,std::string,std::string);
    
    void c_FFVT(std::vector<Single_Vertex>&,int&);
    void c_VVT(std::vector<Single_Vertex>&,int&);
    void c_VVVT(std::vector<Single_Vertex>&,int&);
    
    void c_FFV(std::vector<Single_Vertex>&,int &)  {}
    void c_FFS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &)  {}
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVS(std::vector<Single_Vertex>&,int &)  {}
    void c_VVV(std::vector<Single_Vertex>&,int &)  {}
    void c_VVVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSVV(std::vector<Single_Vertex>&,int &) {}
    void c_SSSS(std::vector<Single_Vertex>&,int &) {}

    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_SSST(std::vector<Single_Vertex>&,int&)  {} 
  };
}

#endif


