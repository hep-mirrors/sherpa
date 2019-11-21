#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE SHERPA::Userhook_Base
#define PARAMETER_TYPE SHERPA::Userhook_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace SHERPA;

Userhook_Base::Userhook_Base(const std::string &name):
  m_name(name)
{
}

Userhook_Base::~Userhook_Base()
{
}

void Userhook_Base::Finish()
{
}

void Userhook_Base::CleanUp()
{
}
