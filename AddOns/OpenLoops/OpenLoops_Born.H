#ifndef AddOns_OpenLoops_OpenLoops_Born_H
#define AddOns_OpenLoops_OpenLoops_Born_H

#include "PHASIC++/Process/Tree_ME2_Base.H"

namespace MODEL {
  class Model_Base;
}

namespace OpenLoops {

  enum AmplitudeType {
    Tree = 1,
    Loop2 = 12
  };

  class OpenLoops_Born: public PHASIC::Tree_ME2_Base {

    int m_ol_id;
    AmplitudeType m_amplitudetype;
    double m_symfac;

  public:

    OpenLoops_Born(const PHASIC::Process_Info &,
                   const ATOOLS::Flavour_Vector &,
                   int ol_id,
                   AmplitudeType);

    ~OpenLoops_Born() {}

    double Calc(const ATOOLS::Vec4D_Vector& momenta);

    int OrderQCD(const int &id=-1);
    int OrderEW(const int &id=-1);

  };

}

#endif
