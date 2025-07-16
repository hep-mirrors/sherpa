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

  class Hadron_Decay_Handler : public Decay_Handler_Base {
    HADRONS::Mixing_Handler* p_mixinghandler;

    long int bs, bs_cs_semilep, bs_ccss, bs_ccds, bs_cuds;
    long int bs_cs_semilep_rej, bs_ccss_rej, bs_ccds_rej, bs_cuds_rej;

  public :
    Hadron_Decay_Handler();
    ~Hadron_Decay_Handler();

    /*!
      \brief Chooses the appropriate decay handler for <var>inpart</var>
      to create the decay blob for inpart. The decay blob is only a stub
      so far without kinematics and without outparticles.
      */
    void CreateDecayBlob(ATOOLS::Particle* part);
    
    void TreatInitialBlob(ATOOLS::Blob* blob,
                          METOOLS::Amplitude2_Tensor* amps,
                          const ATOOLS::Particle_Vector& origparts);
    METOOLS::Amplitude2_Tensor*
    FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);
    METOOLS::Decay_Matrix*
    FillDecayTree(ATOOLS::Blob * blob, METOOLS::Spin_Density* s0);

    /*!
      \brief If a fragmentation blob is encountered, it is checked whether it stems from a
      partonic hadron decay. If it does, and it produced an exclusive decay channel
      of that decaying hadron, it is rejected and redone, as to not spoil that
      hadron's branching ratio.
      */
    bool RejectExclusiveChannelsFromFragmentation(ATOOLS::Blob*);
    /*!
      \brief Calculates the position of the blob in the lab frame according
      to its inparticle's decay time and momentum.
      */
    void SetPosition(ATOOLS::Blob* blob);

    inline double Mass(const ATOOLS::Flavour &fl) const { return fl.HadMass(); }
  };

}

#endif
