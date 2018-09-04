#ifndef CFPSHOWER_Calculators_SF_IF_H
#define CFPSHOWER_Calculators_SF_IF_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_IF : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_IF(const Kernel_Info & info);
    virtual void GeneratePoint(Splitting & split) const {
      split.SetZ();
      split.Setphi();
    }
    virtual bool InitKinematics(Splitting & split) const;
    virtual int  Construct(Splitting & split) const;
  };
}

#endif
