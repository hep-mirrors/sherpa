#ifndef Photon_Remnant_H
#define Photon_Remnant_H

#include "Remnant_Base.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class Photon_Remnant: public Remnant_Base {
  public:
    
    // constructor
    Photon_Remnant(const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
