#ifndef CFPSHOWER_Calculators_Kinematics_FF2_H
#define CFPSHOWER_Calculators_Kinematics_FF2_H

#include "CFPSHOWER++/Calculators/Kinematics_Base.H"

namespace CFPSHOWER {
  class Kinematics_FF2_Soft : public Kinematics_Base {
  protected:
    std::vector<ATOOLS::Vec4D> m_moms;
    ATOOLS::Poincare           m_newsystem;
    ATOOLS::Vec4D              m_momsum, m_pplus, m_pminus;
    double                     m_momscale;
    
    double CalculateY(Splitting & split);
    bool   ConstructSystem(Splitting & split, Configuration & config);
    bool   KinCheck(Splitting & split, Configuration & config);
    bool   CalculateTransverseMomentum(Splitting & split);
    bool   CalculateRescaleFactor();
    void   RescaleAndBoostMomenta(Splitting & split, Configuration & config);
  public:
    Kinematics_FF2_Soft();
    ~Kinematics_FF2_Soft();

    bool operator()(Splitting & split,Configuration & config);
    bool UpdateSystem(Splitting & split,Configuration & config);
    void CalculateWeight(Splitting & split);
  };
}

#endif
