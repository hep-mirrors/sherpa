#ifndef PHASIC__EWSudakov__Comix_Interface_H
#define PHASIC__EWSudakov__Comix_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "COMIX/Main/Single_Process.H"

namespace PHASIC {

  class Comix_Interface {
  public:
    void FillSpinAmplitudes(std::vector<Spin_Amplitudes>&,
                            ATOOLS::Cluster_Amplitude&) const;
  };
} 
#endif
