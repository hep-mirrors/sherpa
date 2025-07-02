#ifndef RECOMBINATION_Main_Recombination_Handler_H
#define RECOMBINATION_Main_Recombination_Handler_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace RECOMBINATION {
  class Recombination_Handler {
  private:
    bool m_on;
 public:
    Recombination_Handler(const bool & on);
    ~Recombination_Handler();
    
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);

    void Initialize();
    void Reset();
  };
}

#endif
