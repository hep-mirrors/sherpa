#ifndef CFPSHOWER_Calculators_SF_FF2_Soft_H
#define CFPSHOWER_Calculators_SF_FF2_Soft_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF2_Soft : public SF_Base {
  protected:
    double                     m_Q2, m_z[2], m_y;
    std::vector<ATOOLS::Vec4D> m_allmoms;
    ATOOLS::Poincare           m_newsystem;
    ATOOLS::Vec4D              m_momsum, m_pplus, m_pminus;
    double                     m_momscale;
    
    double CalculateY(Splitting & split);
    bool   ConstructSystem(Splitting & split, Configuration & config);
    bool   KinCheck(Splitting & split, Configuration & config);
    /*
    bool   CalculateTransverseMomentum(Splitting & split);
    bool   CalculateRescaleFactor();
    void   RescaleAndBoostMomenta(Splitting & split, Configuration & config);
    */
  public:
    SF_FF2_Soft(const Kernel_Info & info) : SF_Base(info) {
      m_logtype = log_type::soft;
      m_type    = kernel_type::FF;
    }
    bool Construct(Splitting & split,Configuration & config);
    bool UpdateSystem(Splitting & split,Configuration & config);
    void CalculateInvariants(Splitting & split);
    void CalculateJacobean(Splitting & split);
  };
}

#endif
