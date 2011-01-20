#ifndef HiddenValley_h
#define HiddenValley_h

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MyComplex.H"

namespace MODEL {
  class Standard_Model;
  class HiddenValley : public Model_Base {
  private :
    Standard_Model *p_sm;
    void ParticleInit();
    void FillSpectrum(PDF::ISR_Handler *const isr);
  public :
    HiddenValley(std::string,std::string,bool);
    ~HiddenValley();
    bool ModelInit(PDF::ISR_Handler *const isr);
  };
}

#endif


