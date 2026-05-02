// HADRON_RESCATTERING/Main/Deuteron_Coalescence.H
#ifndef DEUTERON_COALESCENCE_H
#define DEUTERON_COALESCENCE_H

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"
#include "HADRON_RESCATTERING/XSecs/BaryonBaryon.H"
namespace HADRON_RESCATTERING {

  enum class CoalescenceModel { HardDisk, Exponential, CrossSectionModel };

  class Deuteron_Coalescence {

  private:
    struct FormationResult {
      bool    formed    = false;
      double  P_total   = 0.0;
      std::string  process   = "none";  // "radiative", "single_pion", "double_pion"
      double  sqrtS     = 0.0;
      double  p_rel     = 0.0;
      bool    isAnti    = false;
      ATOOLS::Particle* part1   = nullptr;
      ATOOLS::Particle* part2   = nullptr;
      ATOOLS::Vec4D P_vec;
    };
    double           m_b_coal;
    double           m_p_coal;
    CoalescenceModel m_model;
    BaryonBaryon m_bb; // for cross section model
    ATOOLS::Particle  m_deuteron;        
    FormationResult CheckPair(ATOOLS::Particle* A, ATOOLS::Particle* B);
    FormationResult EvaluateFormation(FormationResult res,const double& b2_fm2);
    ATOOLS::Blob* BuildBlob(const FormationResult& res);

  public:
  
    Deuteron_Coalescence(double b_coal, double p_coal, CoalescenceModel model);
    ATOOLS::Blob* operator()(ATOOLS::Particle* A,
                             ATOOLS::Particle* B,
                             const double& b2_fm2);
  };
}
#endif