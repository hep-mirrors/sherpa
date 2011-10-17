#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class None_Calculator: public Color_Calculator {
  public:

    inline None_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) {}

    std::string Label() const
    {
      return "None";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      m_stat=true;
      return true;
    }

    bool Evaluate(const CObject *a,const CObject *b,
		  const CObject *e)
    {
      m_stat=true;
      return true;
    }

    void AddJ(CObject *const j)
    {
      p_v->AddJ(j);
    }

  };// end of class None_CC

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(None_C_Getter,"None",Color_Calculator,Vertex_Key);

Color_Calculator *None_C_Getter::operator()(const Vertex_Key &key) const
{
  return new None_Calculator(key);
}

void None_C_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"identity";
}

