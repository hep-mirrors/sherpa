#ifndef Interaction_Model_EHC_h
#define Interaction_Model_EHC_h

#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace AMEGIC {  
  class Interaction_Model_EHC : public Interaction_Model_Base {
    Kabbala geff,ghgg,g1,g2,g3,sintW,costW,PL,PR,M_I,root2,vev;
  public:
    Interaction_Model_EHC(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(Single_Vertex *,int &) {}
    void c_FFS(Single_Vertex *,int &) {}
    void c_VVV(Single_Vertex *,int &) {}
    void c_VVS(Single_Vertex *,int &);
    void c_SSV(Single_Vertex *,int &) {}
    void c_SSS(Single_Vertex *,int &) {}
    void c_VVVV(Single_Vertex *,int &);
    void c_SSVV(Single_Vertex *,int &){}
    void c_SSSS(Single_Vertex *,int &){}

    void c_FFT(Single_Vertex  *,int&)  {}
    void c_VVT(Single_Vertex  *,int&)  {} 
    void c_SST(Single_Vertex  *,int&)  {} 
    void c_VVVT(Single_Vertex *,int&)  {} 
    void c_FFVT(Single_Vertex *,int&)  {} 
    void c_SSST(Single_Vertex *,int&)  {} 
  };
}

#endif


