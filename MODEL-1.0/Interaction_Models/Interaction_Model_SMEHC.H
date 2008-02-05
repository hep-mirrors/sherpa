#ifndef Interaction_Model_SMEHC_h
#define Interaction_Model_SMEHC_h

#include "Interaction_Model_SM.H"
#include "Interaction_Model_EHC.H"

namespace MODEL {  
  class Interaction_Model_SMEHC : public Interaction_Model_Base {
    Interaction_Model_SM        * p_mosm;
    Interaction_Model_EHC       * p_mohl;
  public:
    Interaction_Model_SMEHC(MODEL::Model_Base *,std::string,std::string);
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
    ~Interaction_Model_SMEHC();
};
}
#endif
