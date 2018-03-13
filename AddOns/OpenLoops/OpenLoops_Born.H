#ifndef AddOns_OpenLoops_OpenLoops_Born_H
#define AddOns_OpenLoops_OpenLoops_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "AddOns/OpenLoops/OpenLoops_Interface.H"

namespace OpenLoops {

  class OpenLoops_Born: public PHASIC::Tree_ME2_Base {

    int m_ol_id;
    AmplitudeType m_amplitudetype;
    double m_symfac;
    int m_order_ew;
    int m_order_qcd;

  public:

    OpenLoops_Born(const PHASIC::External_ME_Args& args,
                   int ol_id,
                   AmplitudeType);

    ~OpenLoops_Born() {}

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

    int OrderQCD(const int &id=-1) const
      { return m_order_qcd; }

    int OrderEW(const int &id=-1) const
      { return m_order_ew; }

  };

}

#endif
