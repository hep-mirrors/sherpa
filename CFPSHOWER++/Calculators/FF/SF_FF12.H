#ifndef CFPSHOWER_Calculators_SF_FF12_H
#define CFPSHOWER_Calculators_SF_FF12_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF12 : public SF_Base {
  protected:
    double  m_z, m_y, m_ztilde, m_sij;
    virtual double Jacobean(const Splitting & split) const;
  public:
    SF_FF12(const Kernel_Info & info);
    virtual bool Construct(Splitting & split);
  };
}

#endif
