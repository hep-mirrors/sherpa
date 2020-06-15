#ifndef CFPSHOWER_Calculators_SF_FF12_H
#define CFPSHOWER_Calculators_SF_FF12_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF12 : public SF_Base {
  protected:
    virtual double SoftCorrelation(const Splitting * split) const { return 0.; }; 
  public:
    SF_FF12(const Kernel_Info & info);
  };
}

#endif
