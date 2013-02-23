#include "METOOLS/Explicit/Lorentz_Calculator.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class T_Calculator: public Color_Calculator {
  private:

    const CObject *p_a, *p_b;

    int m_type, m_singlet, m_match;

  public:

    inline T_Calculator(const Vertex_Key &key): 
      Color_Calculator(key), m_singlet(false), m_match(true)
    {
      m_cpl=Complex(sqrt(0.5),0.0);
      Flavour fla(key.p_a?key.p_a->Flav():Flavour(kf_gluon));
      Flavour flb(key.p_b?key.p_b->Flav():Flavour(kf_gluon));
      if (key.p_a==NULL && key.p_b==NULL)
	THROW(fatal_error,"Invalid key");
      if (fla.StrongCharge()==3) {
	if (flb.StrongCharge()==-3) m_type=1;
	else m_type=3;
      }
      else if (fla.StrongCharge()==-3) {
	if (flb.StrongCharge()==3) m_type=2;
	else m_type=4;
      }
      else {
	if (flb.StrongCharge()==3) m_type=5;
	else m_type=6;
      }
    }

    std::string Label() const
    {
      return "T";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      p_a=a; 
      p_b=b;
      switch (m_type) {
      case 1:
	m_singlet=(*a)(0)==(*b)(1) && (*a)(0)<=s_cimax;
	m_stat=true;
	return m_stat;
      case 2:
	m_singlet=(*a)(1)==(*b)(0) && (*a)(1)<=s_cimax;
	m_stat=true;
	return m_stat;
      case 3:
	m_singlet=(*b)(0)==(*b)(1) && (*b)(0)>0 && (*b)(0)<=s_cimax;
	m_match=(*a)(0)==(*b)(1);
	m_stat=m_singlet || m_match;
	return m_stat;
      case 4:
	m_singlet=(*b)(0)==(*b)(1) && (*b)(0)>0 && (*b)(0)<=s_cimax;
	m_match=(*a)(1)==(*b)(0);
	m_stat=m_singlet || m_match;
	return m_stat;
      case 5:
	m_singlet=(*a)(0)==(*a)(1) && (*a)(0)>0 && (*a)(0)<=s_cimax;
	m_match=(*b)(0)==(*a)(1);
	m_stat=m_singlet || m_match;
	return m_stat;
      case 6:
	m_singlet=(*a)(0)==(*a)(1) && (*a)(0)>0 && (*a)(0)<=s_cimax;
	m_match=(*b)(1)==(*a)(0);
	m_stat=m_singlet || m_match;
	return m_stat;
      }
      m_stat=false;
      return m_stat;
    }

    void AddJ(CObject *const j)
    {
      switch (m_type) {
      case 1:
	(*j)(0)=(*p_a)(0);
	(*j)(1)=(*p_b)(1);
	break;
      case 2:
	(*j)(0)=(*p_b)(0);
	(*j)(1)=(*p_a)(1);
	break;
      case 3:
	if (m_match) (*j)(0)=(*p_b)(0);
	else (*j)(0)=(*p_a)(0);
	break;
      case 4:
	if (m_match) (*j)(1)=(*p_b)(1);
	else (*j)(1)=(*p_a)(1);
	break;
      case 5:
	if (m_match) (*j)(0)=(*p_a)(0);
	else (*j)(0)=(*p_b)(0);
	break;
      case 6:
	if (m_match) (*j)(1)=(*p_a)(1);
	else (*j)(1)=(*p_b)(1);
	break;
      }
      if (m_singlet) {
	if (m_type<3) {
#ifndef USING__Explicit_OneOverNC_Sum
	  if (p_v->JC()->Out().empty())
#endif
	  {
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
