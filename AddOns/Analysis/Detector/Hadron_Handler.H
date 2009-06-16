#ifndef Analysis_Detector_Hadron_Handler_H
#define Analysis_Detector_Hadron_Handler_H

#include "AddOns/Analysis/Detector/Particle_Smearer_Base.H"
#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Hadron_Handler : public Particle_Smearer_Base {
  private:
    bool   m_fullECal;
    double m_Ecrit, m_exponent, m_ampl, m_kappa, m_lambda, m_threshold_ECal;
    double m_dep_mean_ECal, m_dep_width_ECal;
    double m_dep_mean_HCal, m_dep_width_HCal;
    

    double m_ECalfrac;

    void DetermineTracker();
    void DetermineECal();
    void DetermineHCal();
    void DetermineMuonChambers(){}

    void CalculateEvis(const double E);
    void CalculateDeposits(const double E,const double eta);
  public:
    Hadron_Handler(Primitive_Analysis *);
    ~Hadron_Handler();

    void SetCriticals(const double,const double,const double);
    void SetECalParams(const double,const double);
    void SetHCalParams(const double,const double);
    void SetECalThreshold(const double thr) { m_threshold_ECal = thr; }
    void SetSplitting(const double, const double);
  };
}

#endif
