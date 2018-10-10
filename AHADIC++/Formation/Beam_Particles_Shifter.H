#ifndef AHADIC_Formations_Beam_Particles_Shifter_H
#define AHADIC_Formations_Beam_Particles_Shifter_H

#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Singlet_Tools.H"
#include "AHADIC++/Tools/Constituents.H"
#include <list>

namespace AHADIC {
  class Beam_Particles_Shifter {
  private:
    std::list<Singlet *>      * p_singlets;
    Soft_Cluster_Handler      * p_softclusters;
    std::list<Proto_Particle *> m_beamparts;
    Constituents * p_constituents;
    
    void ExtractBeamParticles();
    bool ShiftBeamParticles();
    bool RescueLightClusters();
  public:
    Beam_Particles_Shifter(std::list<Singlet *> * singlets,
			   Soft_Cluster_Handler * softclusters);
    ~Beam_Particles_Shifter();
    void Init();
    void Reset();

    bool operator()();
  };
}

#endif
