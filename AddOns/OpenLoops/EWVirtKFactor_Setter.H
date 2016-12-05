#ifndef OpenLoops_EWVirtKFactor_H
#define OpenLoops_EWVirtKFactor_H

#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Coupling_Data.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace ATOOLS {
  class NLO_subevt;
}

namespace OpenLoops {

  class EWVirtKFactor_Setter : public PHASIC::KFactor_Setter_Base {
  private:
    ATOOLS::Vec4D_Vector m_p;
    PHASIC::Virtual_ME2_Base* p_ewloop;
    MODEL::Coupling_Map m_cpls;
    double m_deltaew;

    void CopyMomenta();
    void CopyMomenta(const ATOOLS::NLO_subevt& evt);

    void InitEWVirt();
    void CalcEWCorrection();

  public:

    EWVirtKFactor_Setter(const PHASIC::KFactor_Setter_Arguments &args);
    ~EWVirtKFactor_Setter();

    // Default KFactor method
    double KFactor();
    // KFactor for Comix, not yet tested or validated
    double KFactor(const ATOOLS::NLO_subevt& evt);
  };

  class OpenLoops_Interface;
  extern OpenLoops_Interface* s_interface;
}

#endif
