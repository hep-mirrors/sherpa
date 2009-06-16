#ifndef Analysis_Detector_Photon_Handler_H
#define Analysis_Detector_Photon_Handler_H

#include "AddOns/Analysis/Detector/Particle_Smearer_Base.H"
#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Photon_Handler : public Particle_Smearer_Base {
  private:
    void DetermineTracker()     {}
    void DetermineECal();
    void DetermineHCal();
    void DetermineMuonChambers(){}

 public:
    Photon_Handler(Primitive_Analysis *);
    ~Photon_Handler();
  };
}

#endif
