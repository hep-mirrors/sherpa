#ifndef Interaction_Model_Base_h
#define Interaction_Model_Base_h

#include "Model_Base.H"
#include "Vertex.H"
#include "Run_Parameter.H"

namespace AMEGIC {

  struct Single_Vertex; 

  class Interaction_Model_Base {
  protected:  
    MODEL::Model_Base * p_model;
    Vertex            * p_vertex;
    std::string         m_cplscheme,m_yukscheme;
  public:
    Interaction_Model_Base(MODEL::Model_Base *,std::string,std::string);
    Vertex * GetVertex()   { return p_vertex; }
    void     Init_Vertex();

    virtual ~Interaction_Model_Base();

    virtual void c_FFV(Single_Vertex  *,int&) = 0;
    virtual void c_FFS(Single_Vertex  *,int&) = 0;
    virtual void c_VVV(Single_Vertex  *,int&) = 0;
    virtual void c_SSV(Single_Vertex  *,int&) = 0;
    virtual void c_VVS(Single_Vertex  *,int&) = 0;
    virtual void c_SSS(Single_Vertex  *,int&) = 0;
    virtual void c_VVVV(Single_Vertex *,int&) = 0;
    virtual void c_SSVV(Single_Vertex *,int&) = 0;
    virtual void c_SSSS(Single_Vertex *,int&) = 0;
    virtual void c_FFT(Single_Vertex  *,int&) = 0;
    virtual void c_VVT(Single_Vertex  *,int&) = 0;
    virtual void c_SST(Single_Vertex  *,int&) = 0;
    virtual void c_VVVT(Single_Vertex *,int&) = 0;
    virtual void c_FFVT(Single_Vertex *,int&) = 0;
    virtual void c_SSST(Single_Vertex *,int&) = 0;

    virtual int                       ScalarNumber(const std::string);
    virtual double                    ScalarConstant(const std::string);
    virtual ATOOLS::Function_Base * ScalarFunction(const std::string);
    virtual double                    ScalarFunction(const std::string,double);
    virtual ATOOLS::CMatrix         ComplexMatrix(const std::string);
    virtual Complex                   ComplexMatrixElement(const std::string,const int,const int);

    virtual std::string          Name()               { return p_model->Name(); } 
  };
}

extern AMEGIC::Interaction_Model_Base * mo;

#endif







