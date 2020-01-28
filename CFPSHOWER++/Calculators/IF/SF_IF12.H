#ifndef CFPSHOWER_Calculators_SF_IF12_H
#define CFPSHOWER_Calculators_SF_IF12_H

#include "CFPSHOWER++/Calculators/SF_Base.H"

namespace CFPSHOWER {
  class SF_IF12 : public SF_Base {
  protected:
    virtual double Jacobean(const Splitting & split) const;
    virtual double PDFEstimate(const Splitting & split) const;
  public:
    SF_IF12(const Kernel_Info & info);
    virtual bool Construct(Splitting & split,const int & mode=0);
  };
}

#endif
