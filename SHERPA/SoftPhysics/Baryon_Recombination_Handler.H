#ifndef SHERPA_SoftPhysics_Baryon_Recombination_Handler_H
#define SHERPA_SoftPhysics_Baryon_Recombination_Handler_H

#include "RECOMBINATION/Main/Recombination_Handler.H"
#include "ATOOLS/Org/Return_Value.H"

namespace SHERPA {
  class Baryon_Recombination_Handler {
  private:
    bool                                   m_on;
    RECOMBINATION::Recombination_Handler * p_recombination;
  public :
    Baryon_Recombination_Handler();
    ~Baryon_Recombination_Handler();
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const);
    void CleanUp(const size_t & mode=0);
    void Output();
  };
}

#endif
