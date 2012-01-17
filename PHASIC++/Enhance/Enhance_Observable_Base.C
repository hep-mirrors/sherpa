#include "PHASIC++/Enhance/Enhance_Observable_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE PHASIC::Enhance_Observable_Base
#define PARAMETER_TYPE PHASIC::Enhance_Arguments
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace PHASIC;
using namespace ATOOLS;

Enhance_Observable_Base::Enhance_Observable_Base
(const Enhance_Arguments &args): p_proc(args.p_proc)
{
}
