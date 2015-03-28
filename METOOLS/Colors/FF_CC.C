#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Vertex.H"

namespace METOOLS {

  class FF_Calculator: public Color_Calculator {
  private:

    const CObject *p_j[3];

    int  m_mode, m_n[4];
    bool m_maeb, m_meab, m_mbae, m_mbea;

  public:

    inline FF_Calculator(const Vertex_Key &key): 
      Color_Calculator(key) 
    { 
      m_cpl=Complex(0.5,0.0);
      int lid(0), cnt(0), n[4];
      for (int j(0);j<3;++j)
	if (key.p_mv->Color[key.m_n].ParticleArg(j)<0) lid+=j;
	else n[cnt++]=key.p_mv->Color[key.m_n].ParticleArg(j);
      for (int j(0);j<3;++j)
	if (key.p_mv->Color[key.m_n].p_next->ParticleArg(j)<0) lid+=j;
	else n[cnt++]=key.p_mv->Color[key.m_n].p_next->ParticleArg(j);
      if (lid%2) m_cpl=-m_cpl;
      for (int i(0);i<key.p_mv->id.size();++i)
        for (int j(0);j<4;++j)
          if (key.p_mv->id[i]+1==n[j]) m_n[j]=i;
      m_mode=m_n[0]+1==key.p_mv->id.size() ||
        m_n[1]+1==key.p_mv->id.size() ||
        m_n[2]+1==key.p_mv->id.size() ||
        m_n[3]+1==key.p_mv->id.size();
      if (m_mode) {
	if (m_n[0]+1==key.p_mv->id.size() ||
	    m_n[1]+1==key.p_mv->id.size()) {
	  std::swap<int>(m_n[0],m_n[3]);
	  std::swap<int>(m_n[1],m_n[2]);
	}
	if (m_n[2]+1==key.p_mv->id.size()) {
	  std::swap<int>(m_n[2],m_n[3]);
	  m_cpl=-m_cpl;
	}
      }
   }

    std::string Label() const
    {
      return "F*F";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      p_j[0]=j[m_n[0]];
      p_j[1]=j[m_n[1]];
      p_j[2]=j[m_n[2]];
      m_maeb=m_meab=m_mbae=m_mbea=false;
      if ((*p_j[0])(1)==(*p_j[1])(0) &&
	  (*p_j[1])(1)==(*p_j[2])(0)) m_maeb=true;
      if ((*p_j[1])(1)==(*p_j[0])(0) &&
	  (*p_j[0])(1)==(*p_j[2])(0)) m_meab=true;
      if ((*p_j[2])(1)==(*p_j[1])(0) &&
	  (*p_j[1])(1)==(*p_j[0])(0)) m_mbea=true;
      if ((*p_j[2])(1)==(*p_j[0])(0) &&
	  (*p_j[0])(1)==(*p_j[1])(0)) m_mbae=true;
      if (m_maeb && m_meab &&
	  (*p_j[1])(0)==(*p_j[1])(1)) m_maeb=m_meab=false;
      if (m_mbae && m_mbea &&
	  (*p_j[1])(0)==(*p_j[1])(1)) m_mbae=m_mbea=false;
      m_stat=m_maeb || m_meab || m_mbae || m_mbea;
      if (m_mode==0 && m_stat) {
	const CObject *c=j[m_n[3]];
	if (m_maeb) m_maeb=(*p_j[2])(1)==(*c)(0) && (*c)(1)==(*p_j[0])(0);
	if (m_meab) m_meab=(*p_j[2])(1)==(*c)(0) && (*c)(1)==(*p_j[1])(0);
	if (m_mbea) m_mbea=(*p_j[0])(1)==(*c)(0) && (*c)(1)==(*p_j[2])(0);
	if (m_mbae) m_mbae=(*p_j[1])(1)==(*c)(0) && (*c)(1)==(*p_j[2])(0);
	m_stat=m_maeb || m_meab || m_mbae || m_mbea;
	if (m_meab+m_mbae==m_maeb+m_mbea) m_stat=0;
      }
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      if (m_mode==0) {
	if (m_maeb || m_mbea) j->Invert();
	p_v->AddJ(j);
	return;
      }
      if (m_maeb) {
	CObject *c(j->Copy());
	c->Invert();
	(*c)(0)=(*p_j[0])(0);
	(*c)(1)=(*p_j[2])(1);
	p_v->AddJ(c);
      }
      if (m_meab) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[1])(0);
	(*c)(1)=(*p_j[2])(1);
	p_v->AddJ(c);
      }
      if (m_mbae) {
	CObject *c(j->Copy());
	(*c)(0)=(*p_j[2])(0);
	(*c)(1)=(*p_j[1])(1);
	p_v->AddJ(c);
      }
      if (m_mbea) {
	CObject *c(j->Copy());
	c->Invert();
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
