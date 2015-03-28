#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace METOOLS {

  class T_Calculator: public Color_Calculator {
  private:

    const CObject *p_a, *p_b;

    int m_type, m_singlet, m_match, m_n[3];

  public:

    inline T_Calculator(const Vertex_Key &key): 
      Color_Calculator(key), m_singlet(false), m_match(true)
    {
      m_cpl=Complex(sqrt(0.5),0.0);
      int n[3]={key.p_mv->Color[key.m_n].ParticleArg(0),
		key.p_mv->Color[key.m_n].ParticleArg(1),
		key.p_mv->Color[key.m_n].ParticleArg(2)};
      for (int i(0);i<key.p_mv->id.size();++i)
	for (int j(0);j<3;++j)
	  if (key.p_mv->id[i]+1==n[j]) m_n[j]=i;
      m_type=0;
      if (m_n[0]==key.p_mv->id.size()-1) m_type=1;
      if (m_n[2]==key.p_mv->id.size()-1) m_type=2;
      if (m_n[1]==key.p_mv->id.size()-1) m_type=4;
    }

    std::string Label() const
    {
      return "T";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      switch (m_type) {
      case 0: {
	const CObject *g(j[m_n[0]]), *q(j[m_n[1]]), *p(j[m_n[2]]);
	m_singlet=(*g)(0)==(*g)(1) && (*q)(0)==(*p)(1);
	m_match=(*q)(0)==(*g)(1) && (*g)(0)==(*p)(1);
	m_stat=m_singlet || m_match;
	return m_stat;
      }
      case 1:
	p_a=j[m_n[1]]; 
	p_b=j[m_n[2]];
	m_singlet=(*p_a)(0)==(*p_b)(1) && (*p_a)(0)<=s_cimax;
	m_stat=true;
	return m_stat;
      case 2:
	p_a=j[m_n[1]]; 
	p_b=j[m_n[0]];
	m_singlet=(*p_b)(0)==(*p_b)(1) && (*p_b)(0)>0 && (*p_b)(0)<=s_cimax;
	m_match=(*p_a)(0)==(*p_b)(1);
	m_stat=m_singlet || m_match;
	return m_stat;
      case 4:
	p_a=j[m_n[2]]; 
	p_b=j[m_n[0]];
	m_singlet=(*p_b)(0)==(*p_b)(1) && (*p_b)(0)>0 && (*p_b)(0)<=s_cimax;
	m_match=(*p_a)(1)==(*p_b)(0);
	m_stat=m_singlet || m_match;
	return m_stat;
      }
      m_stat=false;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      switch (m_type) {
      case 0:
	(*j)(0)=(*j)(1)=0;
	if (m_singlet) {
	  if (m_match) j->Divide(3.0/2.0);
	  else j->Divide(-3.0);
	}
	p_v->AddJ(j);
	return;
      case 1:
	(*j)(0)=(*p_a)(0);
	(*j)(1)=(*p_b)(1);
	break;
      case 2:
	if (m_match) (*j)(0)=(*p_b)(0);
	else (*j)(0)=(*p_a)(0);
	break;
      case 4:
	if (m_match) (*j)(1)=(*p_b)(1);
	else (*j)(1)=(*p_a)(1);
	break;
      }
      if (m_singlet) {
	if (m_type==1) {
#ifndef USING__Explicit_OneOverNC_Sum
	  if (p_v->JC()->Out().empty())
#endif
	  {
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
	    p_v->AddJ(j);
	    return;
	  }
	}
	else {
	  if (m_match) j->Divide(3.0/2.0);
	  else j->Divide(-3.0);
	}
      }
      p_v->AddJ(j);
    }

  };// end of class T_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(T_Calculator,"T",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,T_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new T_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,T_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"fundamental";
}
