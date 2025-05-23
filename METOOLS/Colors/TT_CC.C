#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class TT_Calculator: public Color_Calculator {
  private:

    const CObject *p_g[2], *p_q[2];

    int m_g[2], m_q[2], m_d, m_n;

  public:

    inline TT_Calculator(const Vertex_Key &key): 
    Color_Calculator(key)
    {
      m_cpl=Complex(0.5,0.0);
      m_n=p_v->V()->id.size()-1;
      if (key.p_mv->id.size()>4) THROW(not_implemented,"Help!");
      int g[2]={key.p_mv->Color[key.m_n].ParticleArg(0),
		key.p_mv->Color[key.m_n].p_next->ParticleArg(0)};
      int q[2]={key.p_mv->Color[key.m_n].ParticleArg(1),
		key.p_mv->Color[key.m_n].p_next->ParticleArg(2)};
      if (q[0]<0 || q[1]<0) {
	std::swap<int>(g[0],g[1]);
	q[1]=key.p_mv->Color[key.m_n].ParticleArg(2);
	q[0]=key.p_mv->Color[key.m_n].p_next->ParticleArg(1);
      }
      if (g[0]<0 || g[1]<0) THROW(fatal_error,"Invalid call");
      for (int i(0);i<key.p_mv->id.size();++i)
        for (int j(0);j<2;++j) {
          if (key.p_mv->id[i]+1==g[j]) m_g[j]=i;
          if (key.p_mv->id[i]+1==q[j]) m_q[j]=i;
	}
      p_g[0]=p_g[1]=p_q[0]=p_q[1]=NULL;
      if (m_g[0]<m_n && m_g[1]<m_n) m_d=m_q[0]<m_n?0:1;
      else m_d=m_g[0]<m_n?0:1;
    }

    std::string Label() const
    {
      return "T*T";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      m_c.clear();
      if (m_g[0]<m_n) p_g[0]=j[m_g[0]];
      if (m_g[1]<m_n) p_g[1]=j[m_g[1]];
      if (m_q[0]<m_n) p_q[0]=j[m_q[0]];
      if (m_q[1]<m_n) p_q[1]=j[m_q[1]];
      if (m_g[0]<m_n && m_g[1]<m_n) {
	if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[m_d])(m_d)==(*p_g[1-m_d])(1-m_d))
	  m_c.push_back(CInfo(m_d?0:(*p_g[1-m_d])(m_d),
			      m_d?(*p_g[1-m_d])(m_d):0,1.0));
	if ((*p_q[m_d])(m_d)==(*p_g[1-m_d])(1-m_d) &&
	    (*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d))
	  m_c.push_back(CInfo(m_d?0:(*p_g[1-m_d])(m_d),
			      m_d?(*p_g[1-m_d])(m_d):0,-1.0/3.0));
	if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[1-m_d])(m_d)==(*p_g[1-m_d])(1-m_d))
	  m_c.push_back(CInfo(m_d?0:(*p_g[m_d])(m_d),
			      m_d?(*p_g[m_d])(m_d):0,-1.0/3.0));
	if ((*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[1-m_d])(m_d)==(*p_g[1-m_d])(1-m_d))
	  m_c.push_back(CInfo(m_d?0:(*p_q[m_d])(m_d),
			      m_d?(*p_q[m_d])(m_d):0,1.0/9.0));
	return m_c.size();
      }
      if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d)) {
	CInfo c(0,0);
	c(1-m_d)=(*p_q[1-m_d])(1-m_d);
	c(m_d)=(*p_g[m_d])(m_d);
	if (c(0)!=c(1)) m_c.push_back(c);
	else {
	  m_c.push_back(CInfo(c(0),c(1),2.0/3.0));
	  for (size_t i(s_cimin);i<=s_cimax;++i)
	    if ((int)i!=c(0)) m_c.push_back(CInfo(i,i,-1.0/3.0));
	}
      }
      if ((*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d)) {
	CInfo c(0,0,-1.0/3.0);
	c(1-m_d)=(*p_q[1-m_d])(1-m_d);
	c(m_d)=(*p_q[m_d])(m_d);
	if (c(0)!=c(1)) m_c.push_back(c);
	else {
	  m_c.push_back(CInfo(c(0),c(1),-2.0/9.0));
	  for (size_t i(s_cimin);i<=s_cimax;++i)
	    if ((int)i!=c(0)) m_c.push_back(CInfo(i,i,1.0/9.0));
	}
      }
      return m_c.size();
    }

  };// end of class TT_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(TT_Calculator,"T*T",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,TT_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new TT_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,TT_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"fundamental*fundamental";
}
