#ifndef CFPSHOWER_Calculators_SF_FI12_H
#define CFPSHOWER_Calculators_SF_FI12_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FI12 : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_FI12(const Kernel_Info & info);
    virtual bool InitKinematics(Splitting & split) const;
    virtual int  Construct(Splitting & split) const;
  };
}

#endif
