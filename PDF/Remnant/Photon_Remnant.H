#ifndef PDF_Remnant_Photon_Remnant_H
#define PDF_Remnant_Photon_Remnant_H

#include "PDF/Remnant/Remnant_Base.H"
#include "PDF/Main/ISR_Handler.H"

namespace PDF {

  class Photon_Remnant: public Remnant_Base {
  public:
    
    // constructor
    Photon_Remnant(const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);
    bool AdjustKinematics();

  };// end of class Remnant_Base

}// end of namespace PDF

#endif
