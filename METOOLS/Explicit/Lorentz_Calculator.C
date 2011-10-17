#include "METOOLS/Explicit/Lorentz_Calculator.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Lorentz_Calculator
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;

Lorentz_Calculator::~Lorentz_Calculator() {}

