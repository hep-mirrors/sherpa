#ifndef SHERPA_SoftPhysics_Lund_Decay_Handler_H
#define SHERPA_SoftPhysics_Lund_Decay_Handler_H

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

namespace SHERPA {

class Lund_Interface;

class Lund_Decay_Handler : public Decay_Handler_Base {

  SHERPA::Lund_Interface* p_lund;

  public :
    Lund_Decay_Handler(SHERPA::Lund_Interface*);
    ~Lund_Decay_Handler();

    /*!
      \brief Chooses the appropriate decay handler for <var>inpart</var>
      to create the decay blob for inpart. The decay blob is only a stub
      so far without kinematics and without outparticles.
      */
    void CreateDecayBlob(ATOOLS::Particle* part);

    METOOLS::Amplitude2_Tensor*
    FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);

    /*!
      \brief Calculates the position of the blob in the lab frame according
      to its inparticle's decay time and momentum.
      */
    void SetPosition(ATOOLS::Blob* blob);

    bool CanDecay(const ATOOLS::Flavour& flav);

    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }
  };

}

#endif
