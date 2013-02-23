#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "ATOOLS/Org/Message.H"

namespace METOOLS {

  class FF_Calculator: public Color_Calculator {
  private:

    const CObject *p_j[3];

    bool m_mab, m_mba, m_sign;

  public:

    inline FF_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    { 
      m_cpl=Complex(0.5,0.0);
    }

    std::string Label() const
    {
      return "F*F";
    }

    bool Evaluate(const CObject *a,const CObject *b,
		  const CObject *e)
    {
      p_j[0]=a;
      p_j[1]=e;
      p_j[2]=b;
      m_mab=m_mba=false;
      if ((*p_j[0])(1)==(*p_j[1])(0) &&
	  (*p_j[2])(0)==(*p_j[1])(1)) m_mab=true;
      if ((*p_j[0])(0)==(*p_j[1])(1) &&
	  (*p_j[2])(1)==(*p_j[1])(0)) m_mba=true;
      if (m_mab && m_mba) // check singlet
	if ((*p_j[0])(1)==(*p_j[2])(0) &&
	    (*p_j[0])(0)==(*p_j[2])(1) &&
	    (*p_j[0])(0)==(*p_j[0])(1)) {
	  m_stat=false;
	  return m_stat;
	}
      m_stat=m_mab || m_mba;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      j->Invert();
      if (m_mab) {
	if (m_mba) {
	  CObject *c(j->Copy());
	  (*c)(0)=(*p_j[2])(0);
	  (*c)(1)=(*p_j[0])(1);
	  p_v->AddJ(c);
	}
	(*j)(0)=(*p_j[0])(0);
	(*j)(1)=(*p_j[2])(1);
      }
      else {
	(*j)(0)=(*p_j[2])(0);
	(*j)(1)=(*p_j[0])(1);
      }
      p_v->AddJ(j);
    }

  };// end of class FF_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(FF_Calculator,"F*F",Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,FF_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new FF_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,FF_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"adjoint";
}
