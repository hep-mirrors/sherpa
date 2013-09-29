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
