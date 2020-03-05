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
  Comix_Interface(Process_Base*, const EWSudakov_Amplitudes&);
  virtual void InitializeProcesses(const Cluster_Amplitude_PM&);
  void FillSpinAmplitudes(std::vector<Spin_Amplitudes>&,
                          ATOOLS::Cluster_Amplitude&) const;

protected:
  Comix_Interface(Process_Base*, const std::string& procname_suffix);
  virtual void InitializeProcess(const Process_Info&);
  const Process_Base* p_proc;
  std::string m_procname_suffix;

private:
  Process_Info CreateProcessInfo(const Cluster_Amplitude*,
                                 const std::string& graph_path);
  void AdaptDifferentialMode();
  NLOTypeStringProcessMap_Map m_apmap;
  int m_differentialmode;

};

} // namespace PHASIC
#endif
