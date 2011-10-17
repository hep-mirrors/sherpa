#ifndef SHERPA_SoftPhysics_Hadron_Decay_Handler_H
#define SHERPA_SoftPhysics_Hadron_Decay_Handler_H

#include <set>
#include <iterator>

#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "SHERPA/Single_Events/Decay_Handler_Base.H"


namespace ATOOLS { 
  class Blob;
  class Blob_List;
}

namespace HADRONS {
  class Mixing_Handler;
}

namespace SHERPA {
  class Soft_Photon_Handler;

  class Hadron_Decay_Handler : public Decay_Handler_Base {
    SHERPA::Soft_Photon_Handler* p_softphotons;
    HADRONS::Mixing_Handler* p_mixinghandler;

  public :
    Hadron_Decay_Handler(std::string path, std::string fragfile);
    ~Hadron_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Particle* part);
    
    void TreatInitialBlob(ATOOLS::Blob* blob,
                          METOOLS::Amplitude2_Tensor* amps,
                          const ATOOLS::Particle_Vector& origparts);
    METOOLS::Amplitude2_Tensor*
    FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);
    METOOLS::Decay_Matrix*
    FillDecayTree(ATOOLS::Blob * blob, METOOLS::Spin_Density* s0);

    bool RejectExclusiveChannelsFromFragmentation(ATOOLS::Blob*);
    void SetPosition(ATOOLS::Blob* blob);
    void AttachExtraQED(ATOOLS::Blob * blob);

    void SetSoftPhotonHandler(Soft_Photon_Handler* sph);
    
    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }
  };

}

#endif
