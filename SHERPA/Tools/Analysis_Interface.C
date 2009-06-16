#include "SHERPA/Tools/Analysis_Interface.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE SHERPA::Analysis_Interface
#define PARAMETER_TYPE SHERPA::Analysis_Arguments
#include "ATOOLS/Org/Getter_Function.C"

using namespace SHERPA;

Analysis_Interface::~Analysis_Interface()
{
}

void Analysis_Interface::CleanUp()
{
}

bool Analysis_Interface::WriteOut()
{
  return true;
}

void Analysis_Interface::ShowSyntax(const int mode)
{
}

