#ifndef AHADIC_Tools_KT_Selector_H
#define AHADIC_Tools_KT_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  class KT_Selector {
  private:
    bool   m_isgluon;
    double m_sigma,m_sigma2,m_sig2,m_offset;
    double WeightFunction(const double & kt);
  public:
    KT_Selector();
    ~KT_Selector();
    void Init(const bool & isgluon);
    
    double operator()(const double & ktmax,const double & M2=-1.);
  };
}

#endif
