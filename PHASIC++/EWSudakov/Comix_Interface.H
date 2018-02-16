#ifndef PHASIC__EWSudakov__Comix_Interface_H
#define PHASIC__EWSudakov__Comix_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "COMIX/Main/Single_Process.H"

namespace PHASIC {

  class EWSudakov_Amplitudes;

  class Comix_Interface {

  public:
    Comix_Interface(Process_Base*, EWSudakov_Amplitudes&);
    void InitializeProcesses(EWSudakov_Amplitudes&);
    void FillSpinAmplitudes(std::vector<Spin_Amplitudes>&,
                            ATOOLS::Cluster_Amplitude&) const;
  private:
    const Process_Base* p_proc;
    NLOTypeStringProcessMap_Map m_apmap;

  };
} 
#endif
