#ifndef PHASIC__EWSudakov__HE_Comix_Interface_H
#define PHASIC__EWSudakov__HE_Comix_Interface_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"

namespace PHASIC {

class HE_Comix_Interface : public Comix_Interface {

public:
  HE_Comix_Interface(Process_Base*, const EWSudakov_Amplitudes&);
  void InitializeProcesses(const Cluster_Amplitude_PM&) override;
  void ResetWithEWParameters(const MODEL::EWParameters&);

protected:
  void InitializeProcess(const Process_Info&) override;

private:
  bool InitializeHighEnergyModel();
  virtual NLOTypeStringProcessMap_Map& ProcessMap() override;
  virtual const NLOTypeStringProcessMap_Map& ProcessMap() const override;
  /// we only want a single model for all interfaces (each process has its own
  /// Sudakov K factor and HE Comix interface), hence we make it static
  static std::unique_ptr<MODEL::Model_Base> p_model_he;
  static NLOTypeStringProcessMap_Map s_apmap_he;
};

} // namespace PHASIC
#endif
