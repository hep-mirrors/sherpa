//bof
//Version: 2 ADICIC++-0.0/2004/09/10

//The main Adicic object.



#ifndef _Adicic_H_
#define _Adicic_H_ _Adicic_H_


#include "Blob_List.H"
#include "Particle_List.H"
#include "Model_Base.H"    //Will give running alpha_S etc.
#include "Cascade_Handler.H"


#define ADICIC_OUTPUT ADICIC_OUTPUT
#undef  ADICIC_OUTPUT





namespace ADICIC {



  class Adicic {
  private:
    std::size_t      m_total;
    std::size_t      m_fail;
    Cascade_Handler* p_handler;
    Cascade          m_cascade;
  public:
    Adicic();
    Adicic(MODEL::Model_Base*);
    ~Adicic();

    void     PrepareCascade();
    Cascade& GetCascade();

    int  PerformShowers();
    bool ExtractPartons(ATOOLS::Blob_List*, ATOOLS::Particle_List* =NULL);
  };



  inline void     Adicic::PrepareCascade() { assert(m_cascade.Clear());}
  inline Cascade& Adicic::GetCascade()     { return m_cascade;}



}    //eo namespace ADICIC





#endif    //eo _Adicic_H_



//eof
