#ifndef Interaction_Model_EW_Grav_h
#define Interaction_Model_EW_Grav_h

#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace AMEGIC {  
  class Interaction_Model_EW_Grav : public Interaction_Model_Base {
    Kabbala    g1,g2,sintW,costW,kap,om,vev;
    Kabbala    PL,PR,M_I,root2,num2,num4,num15;
    Kabbala K_CKM(short int i,short int j);  
    Kabbala conj_K_CKM(short int i,short int j);
 public:
    Interaction_Model_EW_Grav(MODEL::Model_Base *,std::string,std::string);

    void c_FFT(Single_Vertex*,int&);
    void c_VVT(Single_Vertex*,int&);
    void c_SST(Single_Vertex*,int&);
    void c_VVVT(Single_Vertex*,int&);
    void c_FFVT(Single_Vertex*,int&);
    void c_SSST(Single_Vertex*,int&);

    void c_FFV(Single_Vertex *,int &)  {}
    void c_FFS(Single_Vertex *,int &)  {}
    void c_VVV(Single_Vertex *,int &)  {}
    void c_VVS(Single_Vertex *,int &)  {}
    void c_SSV(Single_Vertex *,int &)  {}
    void c_SSS(Single_Vertex *,int &)  {}
    void c_VVVV(Single_Vertex *,int &) {}
    void c_SSVV(Single_Vertex *,int &) {}
    void c_SSSS(Single_Vertex *,int &) {}
  };
}

#endif





