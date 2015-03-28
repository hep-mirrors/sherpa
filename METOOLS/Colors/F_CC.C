#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

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
      for (int i(0);i<key.p_mv->id.size();++i)
	for (int j(0);j<3;++j)
	  if (key.p_mv->id[i]+1==
	      key.p_mv->Color[key.m_n].ParticleArg(j)) m_n[j]=i;
      m_mode=m_n[0]+1==key.p_mv->id.size() ||
	m_n[1]+1==key.p_mv->id.size() ||
	m_n[2]+1==key.p_mv->id.size();
      if (m_mode)
	while (m_n[2]+1!=key.p_mv->id.size())
	  for (int l(m_n[2]),i(2);i>=0;--i) m_n[i]=i?m_n[i-1]:l;
    }

    std::string Label() const
    {
      return "F";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      p_a=j[m_n[0]];
      p_b=j[m_n[1]];
      m_mab=(*p_a)(0)==(*p_b)(1);
      m_mba=(*p_a)(1)==(*p_b)(0);
      m_stat=(m_mab||m_mba)&&
	!((*p_a)(0)==(*p_a)(1) && (*p_b)(0)==(*p_b)(1));
      if (m_mode==0 && m_stat) {
	const CObject *c=j[m_n[2]];
	if (m_mab) m_mab=(*p_b)(0)==(*c)(1) && (*c)(0)==(*p_a)(1);
	if (m_mba) m_mba=(*p_b)(1)==(*c)(0) && (*c)(1)==(*p_a)(0);
	m_stat=m_mab||m_mba;
      }
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      if (m_mode==0) {
	if (m_mab) j->Invert();
      }
      else {
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
