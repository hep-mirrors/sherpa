#ifndef PHASIC__EWSUDAKOV__KFactor__H
#define PHASIC__EWSUDAKOV__KFactor__H

#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/EWSudakov/Clustered_EWSudakov_Calculator.H"

namespace PHASIC {
  class Sudakov_KFactor: public KFactor_Setter_Base {
  private:
    Clustered_EWSudakov_Calculator m_calc;
    ATOOLS::EWSudakov_Log_Corrections_Map m_corrections_map;
    double m_maxweight;
  public:
    Sudakov_KFactor(const KFactor_Setter_Arguments &args);
    double KFactor(const int mode=0);
    double KFactor(const ATOOLS::NLO_subevt &evt);
    const ATOOLS::EWSudakov_Log_Corrections_Map& CorrectionsMap() const
    { return m_corrections_map; }
    const std::set<ATOOLS::EWSudakov_Log_Type>& ActiveLogTypes() const
    { return m_calc.ActiveLogTypes(); }
  };
}

#endif
