#ifndef Four_Particle_Amplitudes_h
#define Four_Particle_Amplitudes_h

#include "Partial_Amplitude_Base.H"
#include "Polarization_Tools.H"

namespace HELICITIES {

  class VSSS : public Partial_Amplitude_Base {
  private:
  public:
    VSSS(ATOOLS::Flavour *fl,int *i,bool *out);
    ~VSSS() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

}

#endif
