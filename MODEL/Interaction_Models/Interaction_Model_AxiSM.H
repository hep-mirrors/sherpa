#ifndef Interaction_Model_AxiSM_h
#define Interaction_Model_AxiSM_h

#include "Interaction_Model_AxiQCD.H"
#include "Interaction_Model_EW.H"
#include "Interaction_Model_Higgs_SM.H"
#include "Interaction_Model_Base.H"
#include "Kabbala.H"

namespace MODEL {  
  class Interaction_Model_AxiSM : public Interaction_Model_Base {
    Interaction_Model_AxiQCD    * p_moqcd;
    Interaction_Model_EW        * p_moew;
    Interaction_Model_Higgs_SM  * p_mosmh;
  public:
    Interaction_Model_AxiSM(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &)  {}
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);
    void c_FFT(std::vector<Single_Vertex>&,int&)  {}
    void c_VVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SST(std::vector<Single_Vertex>&,int&)  {} 
    void c_VVVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_FFVT(std::vector<Single_Vertex>&,int&)  {} 
    void c_SSST(std::vector<Single_Vertex>&,int&)  {} 
    ~Interaction_Model_AxiSM();    
};
}

#endif
