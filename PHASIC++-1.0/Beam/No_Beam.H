#ifndef No_Beam_h
#define No_Beam_h

#include "Beam_Base.H"


namespace BEAM {
  class No_Beam_at_all : public Beam_Base {
  public:
    No_Beam_at_all(APHYTOOLS::Flavour,double);     
    bool   CalculateWeight(double,double);
    double Weight(APHYTOOLS::Flavour);
    bool   On()           { return 0; }
  };
}
#endif
