#ifndef CFPSHOWER_Calculators_SF_FF12_H
#define CFPSHOWER_Calculators_SF_FF12_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF12 : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_FF12(const Kernel_Info & info);
    virtual bool Construct(Splitting & split,const int & mode=0);
  };
}

#endif
