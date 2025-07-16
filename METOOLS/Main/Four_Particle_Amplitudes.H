#ifndef METOOLS_Main_Four_Particle_Amplitudes_H
#define METOOLS_Main_Four_Particle_Amplitudes_H

#include "METOOLS/Main/Partial_Amplitude_Base.H"
#include "METOOLS/Main/Polarization_Tools.H"

namespace METOOLS {

  class VSSS : public Partial_Amplitude_Base {
  private:
  public:
    VSSS(const ATOOLS::Flavour_Vector& flavs,
         const std::vector<int>& i, const std::vector<bool>& out);
    ~VSSS() {}

    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);
  };

}

#endif
