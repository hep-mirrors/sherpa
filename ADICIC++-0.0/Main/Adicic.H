//bof
//Version: 2 ADICIC++-0.0/2004/08/10

//The main Adicic object.



#ifndef _Adicic_H_
#define _Adicic_H_ _Adicic_H_


#include "Blob_List.H"
#include "Particle_List.H"
#include "Model_Base.H"    //Will give running alpha_S etc.
#include "Chain.H"


#define ADICIC_OUTPUT ADICIC_OUTPUT
#undef  ADICIC_OUTPUT





namespace ADICIC {



  class Adicic {
  private:
    Chain* p_chain;
  public:
    Adicic();
    Adicic(MODEL::Model_Base*);
    ~Adicic();

    void   PrepareChain();
    Chain* GetChain();

    int    PerformShowers();
    bool   ExtractPartons(ATOOLS::Blob_List*, ATOOLS::Particle_List* =NULL);
  };



  inline void Adicic::PrepareChain() { p_chain->Clear();}
  inline Chain* Adicic::GetChain() { return p_chain;}



}    //eo namespace ADICIC





#endif    //eo _Adicic_H_



//eof
