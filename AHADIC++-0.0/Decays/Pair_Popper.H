#ifndef Pair_Popper_H
#define Pair_Popper_H

#include "Hadronisation_Parameters.H"
#include "Constituents.H"
#include "Flavour.H"

namespace AHADIC {
  class Pair_Popper {
  private:
    Constituents * p_constituents;
    double         m_totweight;
  public:
    Pair_Popper();
    ~Pair_Popper() {}
    bool Pop(ATOOLS::Flavour &);
  };
}
#endif
