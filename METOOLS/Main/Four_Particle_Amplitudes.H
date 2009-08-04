#ifndef METOOLS_Main_Four_Particle_Amplitudes_H
#define METOOLS_Main_Four_Particle_Amplitudes_H

#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "METOOLS/Main/Polarization_Tools.H"

namespace METOOLS {

  class VSSS : public Partial_Amplitude_Base {
  private:
  public:
    VSSS(ATOOLS::Flavour *fl,int *i,bool *out);
    ~VSSS() {}

    void operator()(const ATOOLS::Vec4D * moms,const bool anti);
  };

}

#endif
