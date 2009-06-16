#ifndef Analysis_Detector_Electron_Handler_H
#define Analysis_Detector_Electron_Handler_H

#include "AddOns/Analysis/Detector/Particle_Smearer_Base.H"
#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Electron_Handler : public Particle_Smearer_Base {
  private:
    void DetermineTracker();
    void DetermineECal();
    void DetermineHCal();
    void DetermineMuonChambers(){}

    double CalculateEFrac(const double);
  public:
    Electron_Handler(Primitive_Analysis *);
    ~Electron_Handler();
  };
}

#endif
