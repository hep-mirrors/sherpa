#ifndef AHADIC_Tools_Z_Selector_H
#define AHADIC_Tools_Z_Selector_H

#include "ATOOLS/Phys/Flavour.H"

namespace AHADIC {
  class Z_Selector {
  private:
    double WeightFunction(const double & z);
  public:
    Z_Selector();
    ~Z_Selector();
    void Init();
    
    double operator()(const double & zmin=1.e-6,const double & zmax=0.999999);
  };
}

#endif
