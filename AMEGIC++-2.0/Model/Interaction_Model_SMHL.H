#ifndef Interaction_Model_SMHL_h
#define Interaction_Model_SMHL_h

#include "Interaction_Model_QCD.H"
#include "Interaction_Model_EW.H"
#include "Interaction_Model_HL.H"

namespace AMEGIC {  
  class Interaction_Model_SMHL : public Interaction_Model_Base {
    Interaction_Model_QCD       * p_moqcd;
    Interaction_Model_EW        * p_moew;
    Interaction_Model_HL        * p_mohl;
  public:
    Interaction_Model_SMHL(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(Single_Vertex  *,int &);
    void c_VVV(Single_Vertex  *,int &);
    void c_VVVV(Single_Vertex *,int &);
    void c_FFS(Single_Vertex  *,int &);
    void c_VVS(Single_Vertex  *,int &);
    void c_SSV(Single_Vertex  *,int &)  {}
    void c_SSS(Single_Vertex  *,int &);
    void c_SSSS(Single_Vertex *,int &);
    void c_SSVV(Single_Vertex *,int &);
    void c_FFT(Single_Vertex  *,int &)  {}
    void c_VVT(Single_Vertex  *,int &)  {}
    void c_SST(Single_Vertex  *,int &)  {}
    void c_VVVT(Single_Vertex *,int &)  {}
    void c_FFVT(Single_Vertex *,int &)  {}
    void c_SSST(Single_Vertex *,int &)  {}
    ~Interaction_Model_SMHL();
};
}
#endif
