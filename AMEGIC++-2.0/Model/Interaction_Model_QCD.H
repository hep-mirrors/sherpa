#ifndef Interaction_Model_QCD_h
#define Interaction_Model_QCD_h

#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace AMEGIC {  
  class Interaction_Model_QCD : public Interaction_Model_Base {
    Kabbala g3,PL,PR,M_I;    
  public:
    Interaction_Model_QCD(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(Single_Vertex *,int &);
    void c_VVV(Single_Vertex *,int &);
    void c_VVVV(Single_Vertex *,int &);

    void c_FFS(Single_Vertex *,int &)  {}
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


