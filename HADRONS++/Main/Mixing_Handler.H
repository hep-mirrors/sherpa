#ifndef HADRONS_Main_Mixing_Handler_H
#define HADRONS_Main_Mixing_Handler_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <map>

namespace ATOOLS {
  class Blob;
  class Particle;
}

namespace HADRONS {
  class Hadron_Decay_Table;
  class Hadron_Decay_Channel;

  /** 
   * Provides the physics needed to deal with mixing phenomena. This includes
   *  - explicit mixing in the event record
   *  - CP rate asymmetries in a decay channel
   */
  class Mixing_Handler {
    private:
    std::map<ATOOLS::Flavour, bool> m_mixing_switch, m_interference_switch;

    /** 
     * Determines the time relevant for mixing of that particle. Usually this
     * is simply the decay time, but in coherent production it is the time
     * difference between the two coherently produced particles.
     * 
     * @param decayer Decaying particle
     * @param checkforpartstatus How to look for coherent sister particle
     * 
     * @return Mixing time.
     */
      double DetermineMixingTime(ATOOLS::Particle* decayer,
                                 bool checkforpartstatus) const;

    public:
    Mixing_Handler(ATOOLS::Scoped_Settings s);
      ~Mixing_Handler();

    /** 
     * Determine whether an explicit mixing event of the decayer happens
     * and if yes, return the created mixing blob.
     * 
     * @param decayer Decaying particle
     * 
     * @return true if mixing happened, false if not.
     */
    bool PerformMixing(ATOOLS::Particle* decayer) const;

    /** 
     * Sets the CP asymmetries for a certain particle's decay table taking
     * into account the interference between mixing and decay of a particle.
     * 
     * @param decayer Decaying particle
     * @param table Hadron_Decay_Table to adjust for these asymmetries
     * 
     * @return true if decay table needed adjustment, false if not
     */
    Hadron_Decay_Channel* Select(ATOOLS::Particle* decayer,
                                 Hadron_Decay_Table& table) const;
  };
}

#endif
