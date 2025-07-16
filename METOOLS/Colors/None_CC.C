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

    bool Evaluate(const CObject_Vector &j)
    {
      m_c.clear();
      m_c.push_back(CInfo(0,0));
      return true;
    }

  };// end of class None_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(None_Calculator,"None",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,None_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new None_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,None_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"identity";
}
