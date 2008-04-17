#ifndef Photon_Handler_h
#define Photon_Handler_h

#include "Particle_Smearer_Base.H"
#include "Particle_List.H"

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
