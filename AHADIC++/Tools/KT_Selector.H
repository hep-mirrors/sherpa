#ifndef AHADIC_Tools_KT_Selector_H
#define AHADIC_Tools_KT_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  class KT_Selector {
  private:
    double m_ktmax,m_sigma,m_sigma2,m_over;
    double WeightFunction(const double & kt);
  public:
    KT_Selector();
    ~KT_Selector();
    void Init();
    
    double operator()(const double & Emax);
  };
}

#endif
