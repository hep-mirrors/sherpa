#ifndef THRew__KFactor__H
#define THRew__KFactor__H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace THRew {
  class THRew_KFactor: public PHASIC::KFactor_Setter_Base {
  public:

    THRew_KFactor(const PHASIC::KFactor_Setter_Arguments& args);
    double KFactor(const int mode = 0) override;
    double KFactor(const ATOOLS::NLO_subevt& evt) override;
    void CalculateAndFillWeightsMap(ATOOLS::Weights_Map&);
    void ResetWeightsMap(ATOOLS::Weights_Map&);

  private:

    void Calculate();
  };
} // namespace THRew

#endif
