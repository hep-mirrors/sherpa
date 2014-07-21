#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class TT_Calculator: public Color_Calculator {
  private:

    const CObject *p_g[2], *p_q[2];

    int m_g[2], m_q[2], m_m, m_d;

  public:

    inline TT_Calculator(const Vertex_Key &key): 
    Color_Calculator(key)
    {
      m_cpl=Complex(0.5,0.0);
      m_g[0]=key.p_mv->Color[key.m_n].ParticleArg(0);
      m_g[1]=key.p_mv->Color[key.m_n].p_next->ParticleArg(0);
      m_q[0]=key.p_mv->Color[key.m_n].ParticleArg(1);
      m_q[1]=key.p_mv->Color[key.m_n].p_next->ParticleArg(2);
      if (m_q[0]==4 || m_q[1]==4) THROW(fatal_error,"Invalid call");
      p_g[0]=p_g[1]=p_q[0]=p_q[1]=NULL;
      if (m_g[0] && m_g[1]) m_d=m_q[0]?0:1;
      else m_d=m_g[0]?0:1;
    }

    std::string Label() const
    {
      return "T*T";
    }

    bool Evaluate(const CObject *a,const CObject *b,const CObject *e)
    {
      m_m=0;
      if (m_g[0]) p_g[0]=m_g[0]==1?a:(m_g[0]==2?b:e);
      if (m_g[1]) p_g[1]=m_g[1]==1?a:(m_g[1]==2?b:e);
      if (m_q[0]) p_q[0]=m_q[0]==1?a:(m_q[0]==2?b:e);
      if (m_q[1]) p_q[1]=m_q[1]==1?a:(m_q[1]==2?b:e);
      if (m_g[0] && m_g[1]) {
	if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[m_d])(m_d)==(*p_g[1-m_d])(1-m_d)) m_m|=1;
	if ((*p_q[m_d])(m_d)==(*p_g[1-m_d])(1-m_d) &&
	    (*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d)) m_m|=2;
	if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[1-m_d])(m_d)==(*p_g[1-m_d])(1-m_d)) m_m|=4;
	if ((*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d) &&
	    (*p_g[1-m_d])(m_d)==(*p_g[1-m_d])(1-m_d)) m_m|=8;
	return m_stat=m_m;
      }
      if ((*p_q[m_d])(m_d)==(*p_g[m_d])(1-m_d)) m_m|=1;
      if ((*p_g[m_d])(m_d)==(*p_g[m_d])(1-m_d)) m_m|=2;
      return m_stat=m_m;
    }

    void AddOctet(CObject *const j)
    {
      if ((*j)(0)!=(*j)(1)) {
	p_v->AddJ(j);
	return;
      }
      CObject *c(j->Copy()), *d(NULL);
      c->Divide(-3.0);
      int cr((*j)(0));
      for (size_t i(s_cimin);i<=s_cimax;++i) {
	if ((int)i==cr) continue;
	(*c)(0)=(*c)(1)=i;
	if (i<s_cimax-(cr==(int)s_cimax)) d=c->Copy();
	p_v->AddJ(c);
	c=d;
      }
      j->Divide(3.0/2.0);
      p_v->AddJ(j);
    }

    void AddJ(CObject *const j)
    {
      if (m_g[0] && m_g[1]) {
	if (m_m&1) {
	  CObject *c(j->Copy());
	  (*c)(m_d)=(*p_g[1-m_d])(m_d);
	  p_v->AddJ(c);
	}
	if (m_m&2) {
	  CObject *c(j->Copy());
	  (*c)(m_d)=(*p_g[1-m_d])(m_d);
	  c->Divide(-3.0);
	  p_v->AddJ(c);
	}
	if (m_m&4) {
	  CObject *c(j->Copy());
	  (*c)(m_d)=(*p_g[m_d])(m_d);
	  c->Divide(-3.0);
	  p_v->AddJ(c);
	}
	if (m_m&8) {
	  CObject *c(j->Copy());
	  (*c)(m_d)=(*p_q[m_d])(m_d);
	  c->Divide(9.0);
	  p_v->AddJ(c);
	}
      }
      else {
	(*j)(1-m_d)=(*p_q[1-m_d])(1-m_d);
	if (m_m&1) {
	  (*j)(m_d)=(*p_g[m_d])(m_d);
	  AddOctet(j->Copy());
	}
	if (m_m&2) {
	  (*j)(m_d)=(*p_q[m_d])(m_d);
	  j->Divide(-3.0);
	  AddOctet(j->Copy());
	}
      }
      j->Delete();
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
