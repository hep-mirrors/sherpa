#include "METOOLS/Explicit/Color_Calculator.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Color_Calculator
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;

size_t Color_Calculator::s_cimin(1);
size_t Color_Calculator::s_cimax(3);

Color_Calculator::~Color_Calculator() {}

bool Color_Calculator::Evaluate(const CObject_Vector &j)
{
  THROW(fatal_error,"Pure virtual method called");
  return false;
}
