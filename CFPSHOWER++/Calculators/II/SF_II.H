#ifndef CFPSHOWER_Calculators_SF_II_H
#define CFPSHOWER_Calculators_SF_II_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_II : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_II(const Kernel_Info & info);
    virtual void GeneratePoint(Splitting & split) const {
      split.SetZ();
      split.Setphi();
    }
    virtual bool InitKinematics(Splitting & split) const;
    virtual int  Construct(Splitting & split) const;
  };
}

#endif
