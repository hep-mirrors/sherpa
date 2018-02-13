#ifndef PHASIC__EWSUDAKOV__KFactor__H
#define PHASIC__EWSUDAKOV__KFactor__H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/EWSudakov/Sudakov.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

namespace PHASIC {
  class Sudakov_KFactor: public KFactor_Setter_Base {
  private:
    //MODEL::One_Running_AlphaS m_runas;
    MODEL::Running_AlphaQED m_runaqed;
    Sudakov m_ews;
    // ATOOLS::Vec4D_Vector m_mom;
  public:
    Sudakov_KFactor(const KFactor_Setter_Arguments &args);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);
  };
}

#endif
