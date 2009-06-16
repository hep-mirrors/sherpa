#ifndef HELICITIES_Main_Four_Particle_Amplitudes_H
#define HELICITIES_Main_Four_Particle_Amplitudes_H

#include "HELICITIES/Main/Partial_Amplitude_Base.H"
#include "HELICITIES/Main/Polarization_Tools.H"

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
