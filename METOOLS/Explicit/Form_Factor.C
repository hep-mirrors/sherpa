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

namespace METOOLS {

  class FFNone: public Form_Factor {
  public:
    FFNone(const Vertex_Key &key): Form_Factor("1",key) {}
    double FF() { return 1; }
  };// end of class FFNone

}// end of namespace METOOLS

DECLARE_GETTER(FFNone,"1",Form_Factor,Vertex_Key);
Form_Factor *Getter<Form_Factor,Vertex_Key,FFNone>::
operator()(const Vertex_Key &args) const
{ return new FFNone(args); }
void ATOOLS::Getter<Form_Factor,Vertex_Key,FFNone>::
PrintInfo(std::ostream &str,const size_t width) const
{ str<<"1"; }
