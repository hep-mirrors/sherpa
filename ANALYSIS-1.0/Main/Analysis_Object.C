#include "Analysis_Object.H"

#include "Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ANALYSIS::Analysis_Object
#define PARAMETER_TYPE ANALYSIS::Argument_Matrix
#include "Getter_Function.C"

using namespace ANALYSIS;
using namespace ATOOLS;

Analysis_Object::Analysis_Object():
  p_ana(NULL), m_obs(false) {}

Analysis_Object::~Analysis_Object()
{
}

void Analysis_Object::Reset()
{
}

void Analysis_Object::EndEvaluation(double scale)
{
}

void Analysis_Object::Output(const std::string &pname)
{
}

Analysis_Object &Analysis_Object::operator+=(const Analysis_Object &obj)
{
  return *this;
}

void Analysis_Object::SetAnalysis(Primitive_Analysis *ana)
{
  p_ana=ana;
}

void Analysis_Object::SetName(const std::string &name)
{
  m_name=name;
}
