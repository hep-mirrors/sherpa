#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class ${color_name}_Calculator: public Color_Calculator {
  private:

    const CObject *p_a, *p_b;

    int m_type, m_n[2];

  public:

    inline ${color_name}_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    {
      int n[2]={key.p_mv->Color[key.m_n].ParticleArg(0),
		key.p_mv->Color[key.m_n].ParticleArg(1)};
      for (size_t i(0);i<key.p_mv->id.size();++i)
	for (int j(0);j<2;++j)
	  if (key.p_mv->id[i]==n[j]-1) m_n[j]=i;
      if (m_n[0]==key.p_mv->id.size()-1)
	std::swap<int>(m_n[0],m_n[1]);
      m_type=m_n[1]==key.p_mv->id.size()-1;
    }

    std::string Label() const
    {
      return "${color_name}";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      p_a=j[m_n[0]];
      if (m_type==0) {
	p_b=j[m_n[1]];
	m_stat=(*p_a)(0)==(*p_b)(1) && (*p_a)(1)==(*p_b)(0);
	if (!m_stat) m_stat=(*p_a)(0)==(*p_a)(1) && (*p_b)(0)==(*p_b)(1);
	return m_stat;
      }
      m_stat=true;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      if (m_type) {
	(*j)(0)=(*p_a)(0);
	(*j)(1)=(*p_a)(1);
	if ((*j)(0)==(*j)(1)) {
	  CObject *c(j->Copy()), *d(NULL);
	  c->Divide(-3.0);
	  int cr((*p_a)(0));
	  for (size_t i(s_cimin);i<=s_cimax;++i) {
	    if ((int)i==cr) continue;
	    (*c)(0)=(*c)(1)=i;
	    if (i<s_cimax-(cr==(int)s_cimax)) d=c->Copy();
	    p_v->AddJ(c);
	    c=d;
	  }
	  j->Divide(3.0/2.0);
	}
      }
      else {
	if ((*p_a)(0)==(*p_a)(1)) {
	  if ((*p_a)(0)!=(*p_b)(1)) j->Divide(-3.0);
	  else j->Divide(3.0/2.0);
	}
      }
      p_v->AddJ(j);
    }

  };

}

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(${color_name}_Calculator,"${color_name}",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,${color_name}_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new ${color_name}_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,${color_name}_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"${color_name}";
}
