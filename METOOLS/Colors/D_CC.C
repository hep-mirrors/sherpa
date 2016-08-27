#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class D_Calculator: public Color_Calculator {
  private:

    int m_type, m_n[2];

  public:

    inline D_Calculator(const Vertex_Key &key): 
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
      return "D";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      m_c.clear();
      const CObject *a(j[m_n[0]]), *b(j[m_n[1]]);
      if (m_type==0) {
	int match((*a)(0)==(*b)(1) && (*a)(1)==(*b)(0));
	int sing((*a)(0)==(*a)(1) && (*b)(0)==(*b)(1));
	if (!match && !sing) return false;
	m_c.push_back(CInfo(0,0,match?(sing?2.0/3.0:1.0):-1.0/3.0));
      }
      else {
	bool sing((*a)(0)==(*a)(1));
	m_c.push_back(CInfo((*a)(0),(*a)(1),sing?2.0/3.0:1.0));
	if (sing)
	  for (size_t i(s_cimin);i<=s_cimax;++i)
	    if ((int)i!=(*a)(0))
	      m_c.push_back(CInfo(i,i,-1.0/3.0));
      }
      return true;
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
