#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

namespace METOOLS {

  class T_Calculator: public Color_Calculator {
  private:

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
      m_c.clear();
      switch (m_type) {
      case 0: {
	const CObject *g(j[m_n[0]]), *q(j[m_n[1]]), *p(j[m_n[2]]);
	m_singlet=(*g)(0)==(*g)(1) && (*q)(0)==(*p)(1);
	m_match=(*q)(0)==(*g)(1) && (*g)(0)==(*p)(1);
	if (!m_singlet && !m_match) return false;
	double c(m_singlet?(m_match?2.0/3.0:-1.0/3.0):1.0);
	m_c.push_back(CInfo(0,0,c));
	return true;
      }
      case 1: {
	CObject *a(j[m_n[1]]), *b(j[m_n[2]]);
	m_singlet=(*a)(0)==(*b)(1) && (*a)(0)<=s_cimax;
#ifndef USING__Explicit_OneOverNC_Sum
	if (p_v->JC()->Out().empty()) 
#endif
	if (m_singlet) {
	  m_c.push_back(CInfo((*a)(0),(*a)(0),2.0/3.0));
	  for (size_t i(s_cimin);i<=s_cimax;++i)
	    if ((int)i!=(*a)(0)) 
	      m_c.push_back(CInfo(i,i,-1.0/3.0));
	  return true;
	}
	m_c.push_back(CInfo((*a)(0),(*b)(1)));
	return true;
      }
      case 2: {
	CObject *a(j[m_n[1]]), *b(j[m_n[0]]);
	m_singlet=(*b)(0)==(*b)(1) && (*b)(0)>0 && (*b)(0)<=s_cimax;
	m_match=(*a)(0)==(*b)(1);
	if (!m_singlet && !m_match) return false;
	double c(m_singlet?(m_match?2.0/3.0:-1.0/3.0):1.0);
	if (m_match) m_c.push_back(CInfo((*b)(0),0,c));
	else m_c.push_back(CInfo((*a)(0),0,c));
	return true;
      }
      case 4: {
	CObject *a(j[m_n[2]]), *b(j[m_n[0]]);
	m_singlet=(*b)(0)==(*b)(1) && (*b)(0)>0 && (*b)(0)<=s_cimax;
	m_match=(*a)(1)==(*b)(0);
	if (!m_singlet && !m_match) return false;
	double c(m_singlet?(m_match?2.0/3.0:-1.0/3.0):1.0);
	if (m_match) m_c.push_back(CInfo(0,(*b)(1),c));
	else m_c.push_back(CInfo(0,(*a)(1),c));
	return true;
      }
      }
      return false;
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
