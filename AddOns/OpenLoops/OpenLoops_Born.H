#ifndef AddOns_OpenLoops_OpenLoops_Born_H
#define AddOns_OpenLoops_OpenLoops_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace OpenLoops {

  class OpenLoops_Born: public PHASIC::Tree_ME2_Base {

    int m_ol_id, m_amptype;
    double m_symfac;

  public:

    OpenLoops_Born(const PHASIC::Process_Info& pi,
                   const ATOOLS::Flavour_Vector& flavs,
                   int ol_id, int amptype);

    ~OpenLoops_Born() {}

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

  };

}

#endif
