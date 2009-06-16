#ifndef MODEL_Interaction_Models_Interaction_Model_Base_H
#define MODEL_Interaction_Models_Interaction_Model_Base_H

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Vertex.H"
#include <vector>

namespace MODEL {

  struct Single_Vertex; 

  struct Interaction_Model_Arguments {
    Model_Base *p_model;
    std::string m_cplscheme, m_yukscheme;
    Interaction_Model_Arguments
    (Model_Base *const model,const std::string &cpl,const std::string &yuk): 
      p_model(model), m_cplscheme(cpl), m_yukscheme(yuk) {}
  };// end of struct Interaction_Model_Arguments

  class Interaction_Model_Base {
  public:

    typedef ATOOLS::Getter_Function
    <Interaction_Model_Base,Interaction_Model_Arguments> 
    Interaction_Model_Getter_Function;

  protected:  
    MODEL::Model_Base * p_model;
    std::string         m_code, m_cplscheme,m_yukscheme;
    bool m_tensors, m_loops, m_agcs;
  public:
    Interaction_Model_Base(const std::string &code,
			   MODEL::Model_Base *,std::string,std::string);

    virtual ~Interaction_Model_Base();

    virtual void c_FFV(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_FFS(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_VVV(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SSV(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_VVS(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SSS(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_VVVV(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SSVV(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SSSS(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_FFT(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_VVT(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SST(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_VVVT(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_FFVT(std::vector<Single_Vertex>&,int&) = 0;
    virtual void c_SSST(std::vector<Single_Vertex>&,int&) = 0;

    virtual int                     ScalarNumber(const std::string);
    virtual double                  ScalarConstant(const std::string);
    virtual Complex                 ComplexConstant(const std::string);
    virtual ATOOLS::Function_Base * ScalarFunction(const std::string);
    virtual double                  ScalarFunction(const std::string,double);
    virtual ATOOLS::CMatrix         ComplexMatrix(const std::string);
    virtual Complex                 ComplexMatrixElement(const std::string,const int,const int);

    virtual std::string          Name();
    inline std::string Code() const { return m_code; }

    inline bool HasLoops() const   { return m_loops;   }
    inline bool HasTensors() const { return m_tensors; }
    inline bool HasAGCs() const    { return m_agcs;    }
  };
}

extern MODEL::Interaction_Model_Base * mo;

#endif







