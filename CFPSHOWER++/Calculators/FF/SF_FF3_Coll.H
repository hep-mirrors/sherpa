#ifndef CFPSHOWER_Calculators_SF_FF3_H
#define CFPSHOWER_Calculators_SF_FF3_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_FF3_Coll : public SF_Base {
  protected:
    double CalculateY(Splitting & split);
    bool   ConstructSystem(Splitting & split);
    bool   KinCheck(Splitting & split);
  public:
    SF_FF3_Coll(const Kernel_Info & info) : SF_Base(info) {
      m_logtype = log_type::coll;
      m_type    = kernel_type::FF;
    }
    bool Construct(Splitting & split,Configuration & config);
    bool UpdateSystem(Splitting & split,Configuration & config);
    void CalculateWeight(Splitting & split);
  };

}

#endif
