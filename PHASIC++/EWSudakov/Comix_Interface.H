#ifndef PHASIC__EWSudakov__Comix_Interface_H
#define PHASIC__EWSudakov__Comix_Interface_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "COMIX/Main/Single_Process.H"
#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"
#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC {

class EWSudakov_Amplitudes;

class Comix_Interface {

public:
  Comix_Interface(Process_Base*, EWSudakov_Amplitudes&);
  void InitializeProcesses(EWSudakov_Amplitudes&);
  void FillSpinAmplitudes(std::vector<Spin_Amplitudes>&,
                          ATOOLS::Cluster_Amplitude&) const;

private:
  Process_Info CreateProcessInfo(const Cluster_Amplitude_UP&,
                                 const std::string& graph_path,
                                 const std::string& suffix="Sudakov");
  void InitializeProcess(const Process_Info&);
  const Process_Base* p_proc;
  NLOTypeStringProcessMap_Map m_apmap, m_apmap_he;

  MODEL::Model_Base          * p_model;
  MODEL::Model_Base          * p_model_he;
  bool InitializeHighEnergyModel();
};
} // namespace PHASIC
#endif
