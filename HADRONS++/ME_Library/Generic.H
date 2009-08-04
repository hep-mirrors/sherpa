#ifndef HADRONS_ME_Library_Generic_H
#define HADRONS_ME_Library_Generic_H
 
#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace METOOLS {
  class Partial_Amplitude_Base;
}

namespace HADRONS {

  class Generic : public HD_ME_Base {
    METOOLS::Partial_Amplitude_Base * p_me;
  public:
    Generic(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name);
    ~Generic();
    void operator()(
                       const ATOOLS::Vec4D  * _p,
                       METOOLS::Spin_Amplitudes * amps);
  };
}
#endif
