#ifndef HADRONS_ME_Library_Generic_H
#define HADRONS_ME_Library_Generic_H
 
#include "HADRONS++/ME_Library/HD_ME_Base.H"

namespace HELICITIES {
  class Partial_Amplitude_Base;
}

namespace HADRONS {

  class Generic : public HD_ME_Base {
    HELICITIES::Partial_Amplitude_Base * p_me;
  public:
    Generic(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name);
    ~Generic();
    void operator()(
                       const ATOOLS::Vec4D  * _p,
                       HELICITIES::Spin_Amplitudes * amps);
  };
}
#endif
