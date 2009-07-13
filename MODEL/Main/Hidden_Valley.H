#ifndef HiddenValley_h
#define HiddenValley_h

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"

namespace MODEL {
  class HiddenValley : public Model_Base {
  private :
    void ParticleInit();
    void FillSpectrum();
  public :
    HiddenValley(std::string,std::string,bool);
    ~HiddenValley();
  };
}

#endif


