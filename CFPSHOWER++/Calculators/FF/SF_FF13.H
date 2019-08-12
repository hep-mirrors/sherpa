#ifndef CFPSHOWER_Calculators_SF_FF13_H
#define CFPSHOWER_Calculators_SF_FF13_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF13 : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_FF13(const Kernel_Info & info);
    virtual bool Construct(Splitting & split);
  };
}

#endif
