#ifndef PHASIC__EWSUDAKOV__KFactor__H
#define PHASIC__EWSUDAKOV__KFactor__H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/EWSudakov/EWSudakov_Calculator.H"

namespace PHASIC {
  class Sudakov_KFactor: public KFactor_Setter_Base {
  private:
    EWSudakov_Calculator m_calc;
    size_t m_check;
  public:
    Sudakov_KFactor(const KFactor_Setter_Arguments &args);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);
  };
}

#endif
