#ifndef MODEL_Interaction_Models_Interaction_Model_FeynRules_h
#define MODEL_Interaction_Models_Interaction_Model_FeynRules_h

#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Math/Kabbala.H"

namespace MODEL {  
  class Interaction_Model_FeynRules : public Interaction_Model_Base {

    ATOOLS::Data_Reader         *p_reader;
    ATOOLS::Algebra_Interpreter *p_algebra;
    std::string m_interactionfile;
  public:
    Interaction_Model_FeynRules(MODEL::Model_Base *,std::string,std::string);
    void c_FFV(std::vector<Single_Vertex>&,int &);
    void c_FFS(std::vector<Single_Vertex>&,int &);
    void c_VVV(std::vector<Single_Vertex>&,int &);
    void c_VVS(std::vector<Single_Vertex>&,int &);
    void c_SSV(std::vector<Single_Vertex>&,int &);
    void c_SSS(std::vector<Single_Vertex>&,int &);
    void c_VVVV(std::vector<Single_Vertex>&,int &);
    void c_SSVV(std::vector<Single_Vertex>&,int &);
    void c_SSSS(std::vector<Single_Vertex>&,int &);
    void c_FFT(std::vector<Single_Vertex>&,int&);
    void c_VVT(std::vector<Single_Vertex>&,int&); 
    void c_SST(std::vector<Single_Vertex>&,int&); 
    void c_VVVT(std::vector<Single_Vertex>&,int&); 
    void c_FFVT(std::vector<Single_Vertex>&,int&); 
    void c_SSST(std::vector<Single_Vertex>&,int&); 
    ~Interaction_Model_FeynRules();    
};
}

#endif
