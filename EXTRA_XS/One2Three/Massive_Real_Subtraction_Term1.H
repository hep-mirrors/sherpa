#pragma once
#include "EXTRA_XS/One2Three/Massive_Real_Subtraction.H"

namespace EXTRAXS {
  class Massive_Real_Subtraction_Term1 : public Massive_Real_Subtraction {
  public:
    Massive_Real_Subtraction_Term1(const std::vector<ATOOLS::Flavour>& flavs,
                               const ATOOLS::Flavour& prop,
                               size_t nonprop, size_t propi, size_t propj)
      : Massive_Real_Subtraction(flavs, prop, nonprop, propi, propj) {}

    void Calculate_real_subtraction(const ATOOLS::Vec4D_Vector& momenta,
                                    bool anti=false) override;
  };
}