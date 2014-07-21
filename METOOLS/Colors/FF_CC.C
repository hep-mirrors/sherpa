#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"
#include "MODEL/Interaction_Models/Color_Function.H"
#include "ATOOLS/Org/Message.H"

namespace METOOLS {

  class FF_Calculator: public Color_Calculator {
  private:

    const CObject *p_j[3];

    int  m_n[4];
    bool m_maeb, m_mbea, m_mabe, m_meba;

  public:

    inline FF_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    { 
      m_cpl=Complex(0.5,0.0);
      int id[4], lp[2], n=0;
      for (int i=0;i<3;++i)
	if ((m_n[n]=key.p_mv->Color[key.m_n].
	     ParticleArg(i)-1)<3) ++n;
	else lp[0]=i;
      for (int i=0;i<3;++i)
	if ((m_n[n]=key.p_mv->Color[key.m_n].
	     p_next->ParticleArg(i)-1)<3) ++n;
	else lp[1]=i;
      for (int nc=0;(nc=m_n[3])>=0;m_n[0]=nc)
	for (int i=4;i>0;--i) m_n[i]=m_n[i-1];
      if ((lp[0]%2)!=(lp[1]%2)) m_cpl=-m_cpl;
   }

    std::string Label() const
    {
      return "F*F";
    }

    bool Evaluate(const CObject *a,const CObject *b,
		  const CObject *e)
    {
      p_j[0]=m_n[0]==0?a:(m_n[0]==1?b:e);
      p_j[1]=m_n[1]==0?a:(m_n[1]==1?b:e);
      p_j[2]=m_n[2]==0?a:(m_n[2]==1?b:e);
      m_maeb=m_mbea=m_mabe=m_meba=false;
      if ((*p_j[0])(1)==(*p_j[1])(0) &&
	  (*p_j[2])(0)==(*p_j[1])(1)) m_maeb=true;
      if ((*p_j[0])(0)==(*p_j[1])(1) &&
	  (*p_j[2])(1)==(*p_j[1])(0)) m_mbea=true;
      if (m_maeb && m_mbea) // check singlet
	if ((*p_j[0])(1)==(*p_j[2])(0) &&
	    (*p_j[0])(0)==(*p_j[2])(1) &&
	    (*p_j[0])(0)==(*p_j[0])(1)) {
	  m_maeb=m_mbea=false;
	}
      if ((*p_j[0])(1)==(*p_j[2])(0) &&
	  (*p_j[1])(0)==(*p_j[2])(1)) m_mabe=true;
      if ((*p_j[0])(0)==(*p_j[2])(1) &&
	  (*p_j[1])(1)==(*p_j[2])(0)) m_meba=true;
      if (m_mabe && m_meba) // check singlet
	if ((*p_j[0])(1)==(*p_j[1])(0) &&
	    (*p_j[0])(0)==(*p_j[1])(1) &&
	    (*p_j[0])(0)==(*p_j[0])(1)) {
	  m_mabe=m_meba=false;
	}
      m_stat=m_maeb || m_mbea || m_mabe || m_meba;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      if (m_mabe) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[0])(0);
	(*c)(1)=(*p_j[1])(1);
	p_v->AddJ(c);
      }
      if (m_meba) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[1])(0);
	(*c)(1)=(*p_j[0])(1);
	p_v->AddJ(c);
      }
      j->Invert();
      if (m_maeb) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[0])(0);
	(*c)(1)=(*p_j[2])(1);
	p_v->AddJ(c);
      }
      if (m_mbea) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[2])(0);
	(*c)(1)=(*p_j[0])(1);
	p_v->AddJ(c);
      }
      j->Delete();
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
