// HADRON_RESCATTERING/Main/Deuteron_Coalescence.H
#ifndef DEUTERON_COALESCENCE_H
#define DEUTERON_COALESCENCE_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"

namespace HADRON_RESCATTERING {

  enum class CoalescenceModel { HardDisk, Exponential };

  class Deuteron_Coalescence {
  private:
    double           m_b_coal;
    double           m_p_coal;
    CoalescenceModel m_model;
    
    ATOOLS::Particle  m_deuteron;        

  public:
  
    Deuteron_Coalescence(double b_coal, double p_coal, CoalescenceModel model);
    ATOOLS::Blob* operator()(ATOOLS::Particle* A,
                             ATOOLS::Particle* B,
                             const double& b2_fm2);
  };
}
#endif