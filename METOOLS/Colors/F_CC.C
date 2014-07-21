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
    int  m_mode, m_n[3];

  public:

    inline F_Calculator(const Vertex_Key &key): 
      Color_Calculator(key), m_mode(0)
    { 
      m_cpl=Complex(0.0,sqrt(0.5));
      m_n[0]=key.m_d?0:key.p_mv->Color[key.m_n].ParticleArg(0);
      m_n[1]=key.m_d?1:key.p_mv->Color[key.m_n].ParticleArg(1);
      m_n[2]=key.m_d?2:key.p_mv->Color[key.m_n].ParticleArg(2);
      m_mode=(m_n[0]>0 && m_n[1]>0 && m_n[2]>0)?1:0;
      if (m_n[1]==0) {
	int n(m_n[0]);
	m_n[0]=m_n[1];
	m_n[1]=m_n[2];
	m_n[2]=n;
      }
      if (m_n[2]==0) {
	int n(m_n[0]);
	m_n[0]=m_n[2];
	m_n[2]=m_n[1];
	m_n[1]=n;
      }
    }

    std::string Label() const
    {
      return "F";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      p_a=m_n[1]==1?a:b;
      p_b=m_n[2]==1?a:b;
      m_mab=(*p_a)(0)==(*p_b)(1);
      m_mba=(*p_a)(1)==(*p_b)(0);
      m_stat=(m_mab||m_mba)&&
	!((*p_a)(0)==(*p_a)(1) && (*p_b)(0)==(*p_b)(1));
      return m_stat;
    }

    bool Evaluate(const CObject *a,const CObject *b,const CObject *e)
    {
      p_a=m_n[1]==1?a:(m_n[1]==2?b:e);
      p_b=m_n[2]==1?a:(m_n[2]==2?b:e);
      m_mab=(*p_a)(0)==(*p_b)(1);
      m_mba=(*p_a)(1)==(*p_b)(0);
      m_stat=(m_mab||m_mba)&&
	!((*p_a)(0)==(*p_a)(1) && (*p_b)(0)==(*p_b)(1));
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

  };// end of class F_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(F_Calculator,"F",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,F_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new F_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,F_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"adjoint";
}
