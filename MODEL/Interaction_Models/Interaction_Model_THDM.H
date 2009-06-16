#ifndef MODEL_Interaction_Models_Interaction_Model_THDM_H
#define MODEL_Interaction_Models_Interaction_Model_THDM_H

#include "MODEL/Interaction_Models/Interaction_Model_Higgs_THDM.H"
#include "MODEL/Interaction_Models/Interaction_Model_SM.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_THDM : public Interaction_Model_Base {
    Interaction_Model_QCD        * p_moqcd;
    Interaction_Model_EW         * p_moew;
    Interaction_Model_Higgs_THDM * p_mohiggs;
  public:
    Interaction_Model_THDM(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &);
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
    ~Interaction_Model_THDM();
  };
}

#endif


