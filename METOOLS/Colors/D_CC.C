#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class D_Calculator: public Color_Calculator {
  private:

    const CObject *p_a, *p_b;

    int m_type;

  public:

    inline D_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    {
      m_type=key.p_a->Flav().Strong()-
	key.p_b->Flav().Strong();
    }

    std::string Label() const
    {
      return "D";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      p_a=a;
      p_b=b;
      if (m_type==0) {
	m_stat=(*a)(0)==(*b)(1) && (*a)(1)==(*b)(0);
	return m_stat;
      }
      m_stat=true;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      switch (m_type) {
      case -1:
	(*j)(0)=(*p_b)(0);
	(*j)(1)=(*p_b)(1);
	break;
      case 1:
	(*j)(0)=(*p_a)(0);
	(*j)(1)=(*p_a)(1);
	break;
      }
      if ((*j)(0) && (*j)(0)==(*j)(1)) {
	CObject *c(j->Copy()), *d(NULL);
	c->Divide(-3.0);
	int cr((*(m_type==1?p_a:p_b))(0));
	for (size_t i(s_cimin);i<=s_cimax;++i) {
	  if ((int)i==cr) continue;
	  (*c)(0)=(*c)(1)=i;
	  if (i<s_cimax-(cr==(int)s_cimax)) d=c->Copy();
	  p_v->AddJ(c);
	  c=d;
	}
	j->Divide(3.0/2.0);
      }
      p_v->AddJ(j);
    }

  };// end of class D_Calculator

  class G_Calculator: public D_Calculator {};

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(D_Calculator,"D",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,D_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new D_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,D_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"delta";
}

DECLARE_GETTER(G_Calculator,"G",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,G_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new D_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,G_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"delta";
}
