#include "PHASIC++/Scales/Core_Scale_Setter.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Shower_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Core_Scale_Setter
#define PARAMETER_TYPE PHASIC::Core_Scale_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

Core_Scale_Setter::Core_Scale_Setter
(const Core_Scale_Arguments &args): p_proc(args.p_proc)
{
}

Core_Scale_Setter::~Core_Scale_Setter()
{
}

namespace PHASIC {

  class Shower_Core_Scale: public Core_Scale_Setter {
  public:

    Shower_Core_Scale(const Core_Scale_Arguments &args):
      Core_Scale_Setter(args) {}

    PDF::CParam Calculate(ATOOLS::Cluster_Amplitude *const ampl)
    {
      return p_proc->Shower()->GetClusterDefinitions()->CoreScale(ampl);
    }

  };// end of class Shower_Core_Scale

}// end of namespace PHASIC

DECLARE_ND_GETTER(Shower_Core_Scale,"SHOWER",
		  Core_Scale_Setter,Core_Scale_Arguments,true);

Core_Scale_Setter *ATOOLS::Getter
<Core_Scale_Setter,Core_Scale_Arguments,Shower_Core_Scale>::
operator()(const Core_Scale_Arguments &args) const
{
  return new Shower_Core_Scale(args);
}

void ATOOLS::Getter<Core_Scale_Setter,Core_Scale_Arguments,Shower_Core_Scale>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"shower core scale"; 
}
