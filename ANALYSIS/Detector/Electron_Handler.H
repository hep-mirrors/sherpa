#ifndef Electron_Handler_h
#define Electron_Handler_h

#include "Particle_Smearer_Base.H"
#include "Particle_List.H"

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
