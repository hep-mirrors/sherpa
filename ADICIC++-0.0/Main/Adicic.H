//bof
//Version: 1 ADICIC++-0.0/2004/07/05

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



  class Adicic_Data_Init {
  private:
    static bool s_status;
  public:
    Adicic_Data_Init() {}
    ~Adicic_Data_Init() {}
  public:
    static inline const bool Status();
    static const bool DoIt();
  };

  inline const bool Adicic_Data_Init::Status() {    //Static.
    return s_status;
  }





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
