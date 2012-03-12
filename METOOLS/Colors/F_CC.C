#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "ATOOLS/Org/Message.H"

namespace METOOLS {

  class F_Calculator: public Color_Calculator {
  private:

    const CObject *p_a, *p_b;

    bool m_mab, m_mba;

  public:

    inline F_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    { 
      m_cpl=Complex(0.0,sqrt(0.5));
    }

    std::string Label() const
    {
      return "F";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      p_a=a;
      p_b=b;
      m_mab=(*a)(0)==(*b)(1);
      m_mba=(*a)(1)==(*b)(0);
      m_stat=(m_mab||m_mba)&&
	!((*a)(0)==(*a)(1) && (*b)(0)==(*b)(1));
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      if (m_mab) {
	if (m_mba) {
	  CObject *c(j->Copy());
	  (*c)(0)=(*p_a)(0);
	  (*c)(1)=(*p_b)(1);
	  p_v->AddJ(c);
	}
	(*j)(0)=(*p_b)(0);
	(*j)(1)=(*p_a)(1);
	j->Invert();
      }
      else {
	(*j)(0)=(*p_a)(0);
	(*j)(1)=(*p_b)(1);
      }
      p_v->AddJ(j);
    }

  };// end of class Color_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(F_C_Getter,"F",Color_Calculator,Vertex_Key);

Color_Calculator *F_C_Getter::operator()(const Vertex_Key &key) const
{
  return new F_Calculator(key);
}

void F_C_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"adjoint";
}

