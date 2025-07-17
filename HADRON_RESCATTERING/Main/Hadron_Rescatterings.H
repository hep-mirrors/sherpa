#ifndef HADRON_RESCATTERING_Main_Hadron_Rescatterings_H
#define HADRON_RESCATTERING_Main_Hadron_Rescatterings_H

#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Org/Return_Value.H"

namespace HADRON_RESCATTERING {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Hadron_Rescatterings {
  private:
    bool                m_on;
    long unsigned int   m_nfails;

    BaryonBaryon *      p_BB;
 public:
    Hadron_Rescatterings(const bool & on);
    ~Hadron_Rescatterings();


    
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);

    void Initialize();
    void Reset();
  };
}

#endif
