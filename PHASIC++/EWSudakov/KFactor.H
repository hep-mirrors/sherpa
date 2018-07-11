#ifndef PHASIC__EWSUDAKOV__KFactor__H
#define PHASIC__EWSUDAKOV__KFactor__H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/EWSudakov/Sudakov.H"

namespace PHASIC {
  class Sudakov_KFactor: public KFactor_Setter_Base {
  private:
    Sudakov m_ews;
    size_t m_check;
  public:
    Sudakov_KFactor(const KFactor_Setter_Arguments &args);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);
  };
}

#endif
