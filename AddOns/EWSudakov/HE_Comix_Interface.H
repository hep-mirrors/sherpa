#ifndef PHASIC__EWSudakov__HE_Comix_Interface_H
#define PHASIC__EWSudakov__HE_Comix_Interface_H

#include "AddOns/EWSudakov/Comix_Interface.H"

namespace MODEL {
  class EWParameters;
}

namespace EWSud {

class HE_Comix_Interface : public Comix_Interface {

public:
  HE_Comix_Interface(PHASIC::Process_Base*, const EWSudakov_Amplitudes&);
  void ResetWithEWParameters(const MODEL::EWParameters&);

protected:
  void InitializeProcesses(const Cluster_Amplitude_PM&);

private:
  bool InitializeHighEnergyModel();
  virtual PHASIC::NLOTypeStringProcessMap_Map& ProcessMap() override;
  virtual const PHASIC::NLOTypeStringProcessMap_Map& ProcessMap() const override;
  /// we only want a single model for all interfaces (each process has its own
  /// EWSud K factor and HE Comix interface), hence we make it static
  static std::unique_ptr<MODEL::Model_Base> p_model_he;
  static PHASIC::NLOTypeStringProcessMap_Map s_apmap_he;
};

} // namespace PHASIC
#endif
