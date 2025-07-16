#include "SHERPA/Tools/Event_Reader_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE SHERPA::Event_Reader_Base
#define PARAMETER_TYPE SHERPA::Input_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace SHERPA;

Event_Reader_Base::~Event_Reader_Base()
{
}
