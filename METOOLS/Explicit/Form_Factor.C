#include "METOOLS/Explicit/Form_Factor.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Form_Factor
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;
using namespace ATOOLS;

Form_Factor::Form_Factor(const std::string &id,const Vertex_Key &key):
  m_id(id), p_v(key.p_v) {}

Form_Factor::~Form_Factor()
{
}

std::ostream &METOOLS::operator<<(std::ostream &str,const Form_Factor &c)
{
  return str<<c.ID();
}
