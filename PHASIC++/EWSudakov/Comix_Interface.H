#ifndef PHASIC__EWSudakov__Comix_Interface_H
#define PHASIC__EWSudakov__Comix_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "COMIX/Main/Single_Process.H"

namespace PHASIC {

  class Comix_Interface {
  public:
    Comix_Interface(Process_Base&, ATOOLS::Cluster_Amplitude&);
    ~Comix_Interface();
    void FillSpinAmplitudes(std::vector<Spin_Amplitudes>&);
  private:
    Process_Base& m_proc;
    ATOOLS::Cluster_Amplitude& m_ampl;
    NLOTypeStringProcessMap_Map m_apmap;
  };
} 
#endif
