#ifndef EWSudakov__Comix_Interface_H
#define EWSudakov__Comix_Interface_H

#include "METOOLS/Main/Spin_Structure.H"
#include "AddOns/EWSud/Amplitudes.H"
#include "PHASIC++/Process/Process_Base.H"

namespace ATOOLS {
class Cluster_Amplitude;
}

namespace EWSud {

class Amplitudes;

class Comix_Interface {

public:
  Comix_Interface(PHASIC::Process_Base*, const Amplitudes&);

  /// does not assume that a flavour-sorted Cluster Amplitude is passed,
  /// spin amplitudes for all spins are cached, and can be reset using
  /// ResetSpinAmplitudeCache
  Complex GetSpinAmplitude(ATOOLS::Cluster_Amplitude&,
                           const std::vector<int>& spins) const;
  void ResetSpinAmplitudeCache() const { m_spinampls.clear(); }

  /// assumes that a flavour-sorted Cluster Amplitude is passed
  void FillSpinAmplitudes(std::vector<METOOLS::Spin_Amplitudes>&,
                          const ATOOLS::Cluster_Amplitude&) const;

  // facilities to generate a process from a given amplitude (via an
  // intermediate Process_Info)
  PHASIC::Process_Base* GetProcess(const ATOOLS::Cluster_Amplitude& ampl) const;
  PHASIC::Process_Info CreateProcessInfo(const ATOOLS::Cluster_Amplitude*,
                                         const std::string& graph_path="");
  PHASIC::Process_Base* InitializeProcess(const PHASIC::Process_Info&);

protected:
  Comix_Interface(PHASIC::Process_Base*, const std::string& procname_suffix);
  void InitializeProcesses(const Cluster_Amplitude_PM&);
  const PHASIC::Process_Base* p_proc;
  std::string m_procname_suffix;

private:
  virtual PHASIC::NLOTypeStringProcessMap_Map& ProcessMap()
  { return s_apmap; }
  virtual const PHASIC::NLOTypeStringProcessMap_Map& ProcessMap() const
  { return s_apmap; }
  void AdaptToProcessColorScheme();
  void SetScales(ATOOLS::Cluster_Amplitude&) const;

  int m_differentialmode;
  static PHASIC::NLOTypeStringProcessMap_Map s_apmap;

  /// caches for spin amplitudes and for permutations and permutation signs due
  /// to re-ordering fermions in the final state
  template<class T>
  using AmplMap = std::map<const ATOOLS::Cluster_Amplitude*, T>;
  mutable AmplMap<std::vector<METOOLS::Spin_Amplitudes>> m_spinampls;
  mutable AmplMap<std::vector<size_t>> m_permutations;
  mutable AmplMap<double> m_permutation_signs;

  static double CalcPermutationSign(std::vector<size_t> perm,
                                    const ATOOLS::Cluster_Amplitude&);
};

} // namespace PHASIC
#endif
