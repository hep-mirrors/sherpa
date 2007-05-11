#ifndef Muon_Handler_h
#define Muon_Handler_h

#include "Particle_Smearer_Base.H"
#include "Particle_List.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Muon_Handler : public Particle_Smearer_Base {
  private:
    double m_dep_mean_ECal, m_dep_width_ECal;
    double m_dep_mean_HCal, m_dep_width_HCal;
    
    void DetermineTracker();
    void DetermineECal();
    void DetermineHCal();
    void DetermineMuonChambers();

  public:
    Muon_Handler(Primitive_Analysis *);
    ~Muon_Handler();

    void SetECalParams(const double,const double);
    void SetHCalParams(const double,const double);
  };
}

#endif
