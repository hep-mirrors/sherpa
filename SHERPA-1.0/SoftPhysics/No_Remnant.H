#ifndef No_Remnant_H
#define No_Remnant_H

#include "Remnant_Base.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class No_Remnant: public Remnant_Base {
  public:
    
    // constructor
    No_Remnant(const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);
    bool AdjustKinematics();

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
