#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ANALYSIS::Analysis_Object
#define PARAMETER_TYPE ANALYSIS::Argument_Matrix
#include "ATOOLS/Org/Getter_Function.C"

using namespace ANALYSIS;
using namespace ATOOLS;

Analysis_Object::Analysis_Object():
  p_ana(NULL), m_isobs(false), m_isdet(false) {}

Analysis_Object::~Analysis_Object() {}

void Analysis_Object::EvaluateNLOcontrib(double value, double ncount)
{
  msg_Error()<<"ERROR virtual function Analysis_Object::EvaluateNLOcontrib called "<<m_name<<std::endl
     <<" not NLO-ready!!"<<m_name<<std::endl;
}
 
void Analysis_Object::EvaluateNLOevt()
{
  msg_Error()<<"ERROR virtual function Analysis_Object::EvaluateNLOevt called "<<m_name<<std::endl
     <<" not NLO-ready!!"<<m_name<<std::endl;
}

void Analysis_Object::Reset() {}

void Analysis_Object::Restore(double scale) {}

void Analysis_Object::EndEvaluation(double scale) {}

void Analysis_Object::Test(const int mode) {}

void Analysis_Object::Output(const std::string &pname) {}

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
