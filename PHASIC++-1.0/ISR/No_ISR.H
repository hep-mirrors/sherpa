#ifndef No_ISR_h
#define No_ISR_h

#include "ISR_Base.H"

namespace ISR {
  class No_ISR_at_all : public ISR_Base {
  public:
    No_ISR_at_all(APHYTOOLS::Flavour);     
    bool   CalculateWeight(double,double);
    double Weight(APHYTOOLS::Flavour);
    bool   On()           { return 0; }
    PDF::PDF_Base * PDF() { return 0; }
  };
}
#endif



