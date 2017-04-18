#ifndef AHADIC_Tools_Z_Selector_H
#define AHADIC_Tools_Z_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  class Z_Selector {
  private:
    double m_alpha[3], m_beta[3], m_gamma[3];
    size_t m_flcnt;
    
    ATOOLS::Flavour m_flav;
    double WeightFunction(const double & z);
    double MassFunction(const double & m2, const double & z);
  public:
    Z_Selector();
    ~Z_Selector();
    void Init();
    
    double operator()(const double & zmin=1.e-6,const double & zmax=0.999999,
		      const ATOOLS::Flavour & flav1=ATOOLS::Flavour(kf_none));
  };
}

#endif
