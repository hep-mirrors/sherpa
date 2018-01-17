#ifndef PHASIC__EWSudakov__Sudakov_H
#define PHASIC__EWSudakov__Sudakov_H

#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"

namespace PHASIC {
  class Comix_Interface;
  class Sudakov {
  private:
    Process_Base *p_proc;
    ATOOLS::Cluster_Amplitude *p_ampl;
    Comix_Interface *p_ci;
  public:
    Sudakov(Process_Base *const p_proc);
    ~Sudakov();
    double EWSudakov(ATOOLS::Vec4D_Vector &mom);
  };
}
#endif
