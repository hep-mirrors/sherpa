#ifndef PHASIC__EWSudakov__HE_Comix_Interface_H
#define PHASIC__EWSudakov__HE_Comix_Interface_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"

namespace PHASIC {

class HE_Comix_Interface : public Comix_Interface {

public:
  HE_Comix_Interface(Process_Base*, EWSudakov_Amplitudes&);
  void InitializeProcesses(EWSudakov_Amplitudes&) override;
  void ResetWithEWParameters(const MODEL::EWParameters&);

private:
  void InitializeProcess(const Process_Info&) override;
  bool InitializeHighEnergyModel();

  /// we only want a single model for all interfaces (each process has its own
  /// Sudakov K factor and HE Comix interface), hence we make it static
  static std::unique_ptr<MODEL::Model_Base> p_model_he;
};

} // namespace PHASIC
#endif
