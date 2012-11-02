#ifndef PHASIC_Scales_Core_Scale_Setter_H
#define PHASIC_Scales_Core_Scale_Setter_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PDF/Main/Cluster_Definitions_Base.H"

namespace PHASIC {

  class Process_Base;

  struct Core_Scale_Arguments {
    Process_Base *p_proc;
    std::string m_scale;
    Core_Scale_Arguments(Process_Base *const proc,
			 const std::string &scale=""):
      p_proc(proc), m_scale(scale) {}
  };// end of class Core_Scale_Arguments

  class Core_Scale_Setter {
  protected:

    Process_Base *p_proc;

  public:

    Core_Scale_Setter(const Core_Scale_Arguments &args);

    virtual ~Core_Scale_Setter();

    virtual PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl) = 0;

  };// end of class Core_Scale_Setter

  typedef ATOOLS::Getter_Function
  <Core_Scale_Setter,Core_Scale_Arguments> Core_Scale_Getter;

}// end of namespace PHASIC

#endif
