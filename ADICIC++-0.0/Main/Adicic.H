#ifndef Adicic_h
#define Adicic_h

#include "Blob_List.H"
#include "Particle_List.H"
#include "Chain.H"
//#include "Model_Base.H" // Will give Running Alpha_S etc.

namespace ADICIC {
  class Adicic {
  private :
    Chain * p_chain;
  public :
    Adicic(); 
    ~Adicic(); 
    int  PerformShowers();
    bool ExtractPartons(ATOOLS::Blob_List *,ATOOLS::Particle_List * = NULL);
    void PrepareChain() { p_chain->Clear(); }

    Chain * GetChain();
  };
  inline Chain * Adicic::GetChain() { return p_chain; }
}
#endif







