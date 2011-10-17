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
    Generic(const ATOOLS::Flavour_Vector& flavs,
            const std::vector<int>& decayindices,
            const std::string& name);
    ~Generic();
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };
}
#endif
