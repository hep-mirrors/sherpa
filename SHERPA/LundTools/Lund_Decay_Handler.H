#ifndef SHERPA_SoftPhysics_Lund_Decay_Handler_H
#define SHERPA_SoftPhysics_Lund_Decay_Handler_H

#include <set>
#include <iterator>

#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "PHASIC++/Decays/Decay_Handler.H"


namespace ATOOLS {
  class Blob;
  class Blob_List;
}

namespace SHERPA {

  class Lund_Interface;

  class Lund_Decay_Handler : public PHASIC::Decay_Handler {

    SHERPA::Lund_Interface* p_lund;

  public :
    Lund_Decay_Handler(SHERPA::Lund_Interface*);
    ~Lund_Decay_Handler();

    void CreateDecayBlob(ATOOLS::Blob_List* bloblist, ATOOLS::Particle* part);
    
    METOOLS::Amplitude2_Tensor*
    FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);

    void SetPosition(ATOOLS::Blob* blob);

    bool CanDecay(const ATOOLS::Flavour& flav);

    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }
    bool DiceMass(ATOOLS::Particle* part, double max);
    inline std::string Name() { return "LundDecays"; }
  };

}

#endif
