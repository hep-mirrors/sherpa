#ifndef REMNANTS_Main_Photon_Remnant_H
#define REMNANTS_Main_Photon_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"

namespace REMNANTS {
  class Photon_Remnant: public Remnant_Base {
  public:
    
    // constructor
    Photon_Remnant(const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);
    bool AdjustKinematics();

  };// end of class Remnant_Base

}// end of namespace REMNANTS

#endif
