#include "SHERPA/Tools/Output_Base.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE SHERPA::Output_Base
#define PARAMETER_TYPE SHERPA::Output_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace SHERPA;

Output_Base::Output_Base(const std::string &name):
  m_name(name), p_eventhandler(NULL)
{
}

Output_Base::~Output_Base()
{
}

void Output_Base::Header()
{
}

void Output_Base::Footer()
{
}

void Output_Base::ChangeFile()
{
}

void Output_Base::SetXS(const double&, const double&)
{
}
