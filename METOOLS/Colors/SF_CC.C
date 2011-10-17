#include "METOOLS/Explicit/Dipole_Color.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class SF_Calculator: public Dipole_Color {
  private:

    int m_ti, m_tk;

  public:

    inline SF_Calculator(const Vertex_Key &key): 
      Dipole_Color(key)
    {
      m_cpl=p_cc->Coupling();
      if (key.p_c->Flav().StrongCharge()!=8)
	THROW(fatal_error,"Invalid call");
      m_ti=key.p_a->Flav().StrongCharge();
      m_tk=key.p_k->Flav().StrongCharge();
      if (m_ti!=8) m_cpl/=sqrt(6.0);
    }

    std::string Label() const
    {
      return "S-F";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      m_ci.clear();
      m_cjk.clear();
      m_stat=true;
      switch (m_ti) {
      case 8:
	switch (m_tk) {
	case -3:
	  if ((*a)(0)!=(*b)(1)) {
	    m_ci.push_back(CInfo((*a)(0),(*b)(1),1,0));
	    m_cjk.push_back(CInfo(0,(*a)(1),1,0));
	  }
	  else {
	    bool s((*a)(0)==(*a)(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(0))) {
		m_ci.push_back(CInfo(i,(*a)(1),0,0));
		m_cjk.push_back(CInfo(0,i,1,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*a)(0),(*b)(1),1,0));
	  }
	  break;
	case 3:
	  if ((*a)(1)!=(*b)(0)) {
	    m_ci.push_back(CInfo((*b)(0),(*a)(1),0,0));
	    m_cjk.push_back(CInfo((*a)(0),0,0,0));
	  }
	  else {
	    bool s((*a)(1)==(*a)(0));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(1))) {
		m_ci.push_back(CInfo((*a)(0),i,1,0));
		m_cjk.push_back(CInfo(i,0,0,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*b)(0),(*a)(1),0,0));
	  }
	  break;
	case 8:
	  if ((*a)(1)!=(*b)(0)) {
	    m_ci.push_back(CInfo((*b)(0),(*a)(1),0,0));
	    m_cjk.push_back(CInfo((*a)(0),(*b)(1),0,0));
	  }
	  else {
	    bool s((*a)(0)==(*a)(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(0))) {
		m_ci.push_back(CInfo((*a)(0),i,1,0));
		m_cjk.push_back(CInfo(i,(*b)(1),0,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*b)(0),(*a)(1),0,0));
	  }
	  if ((*a)(0)!=(*b)(1)) {
	    m_ci.push_back(CInfo((*a)(0),(*b)(1),1,1));
	    m_cjk.push_back(CInfo((*b)(0),(*a)(1),1,1));
	  }
	  else {
	    bool s((*a)(1)==(*a)(0));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(1))) {
		m_ci.push_back(CInfo(i,(*a)(1),0,1));
		m_cjk.push_back(CInfo((*b)(0),i,1,1));
	      }
	    if (!s) m_ci.push_back(CInfo((*a)(0),(*b)(1),1,1));
	  }
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	break;
      default: THROW(fatal_error,"Invalid call");
      }
      return m_stat;
    }

    bool Evaluate(const CObject *a,const CObject *b,const CObject *c)
    {
      m_ci.clear();
      m_cjk.clear();
      m_stat=p_cc->Evaluate(a,b);
      if (!m_stat) return false;
      switch (m_ti) {
      case 8:
	m_ci.push_back(CInfo((*a)(0),(*a)(1),1,0));
	m_ci.push_back(CInfo((*b)(0),(*b)(1),0,1));
	switch (m_tk) {
	case 3:
	  if ((*b)(1)==(*c)(0)) m_cjk.push_back(CInfo((*b)(0),0,0,0));
	  if ((*b)(1)==(*b)(0)) m_cjk.push_back(CInfo((*c)(0),0,0,0,-3.0));
	  if ((*a)(1)==(*c)(0)) m_cjk.push_back(CInfo((*a)(0),0,0,1));
	  if ((*a)(1)==(*a)(0)) m_cjk.push_back(CInfo((*c)(0),0,0,1,-3.0));
	  break;
	case -3:
	  if ((*b)(0)==(*c)(1)) m_cjk.push_back(CInfo(0,(*b)(1),1,0));
	  if ((*b)(0)==(*b)(1)) m_cjk.push_back(CInfo(0,(*c)(1),1,0,-3.0));
	  if ((*a)(0)==(*c)(1)) m_cjk.push_back(CInfo(0,(*a)(1),1,1));
	  if ((*a)(0)==(*a)(1)) m_cjk.push_back(CInfo(0,(*c)(1),1,1,-3.0));
	  break;
	case 8:
	  if ((*b)(0)==(*c)(1)) m_cjk.push_back(CInfo((*c)(0),(*b)(1),1,0));
	  if ((*c)(0)==(*b)(1)) m_cjk.push_back(CInfo((*b)(0),(*c)(1),0,0));
	  if ((*a)(0)==(*c)(1)) m_cjk.push_back(CInfo((*c)(0),(*a)(1),1,1));
	  if ((*c)(0)==(*a)(1)) m_cjk.push_back(CInfo((*a)(0),(*c)(1),0,1));
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	break;
      case -3: {
	switch (m_tk) {
	case 3: {
	  if ((*a)(1)!=(*c)(0)) {
	    m_ci.push_back(CInfo((*c)(0),(*a)(1),0,0));
	    m_cjk.push_back(CInfo((*b)(0),0,0,0));
	  }
	  else {
	    bool s((*b)(0)==(*c)(0));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*b)(0))) {
		m_ci.push_back(CInfo((*b)(0),i,1,0));
		m_cjk.push_back(CInfo(i,0,0,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*c)(0),(*a)(1),0,0));
	  }
	  break;
	}
	case -3: {
	  if ((*b)(0)!=(*c)(1)) {
	    m_ci.push_back(CInfo((*b)(0),(*c)(1),1,0));
	    m_cjk.push_back(CInfo(0,(*a)(1),1,0));
	  }
	  else {
	    bool s((*a)(1)==(*c)(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(1))) {
		m_ci.push_back(CInfo(i,(*a)(1),0,0));
		m_cjk.push_back(CInfo(0,i,1,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*b)(0),(*c)(1),1,0));
	  }
	  break;
	}
	case 8: {
	  if ((*a)(1)!=(*c)(0)) {
	    m_ci.push_back(CInfo((*c)(0),(*a)(1),0,0));
	    m_cjk.push_back(CInfo((*b)(0),(*c)(1),0,0));
	  }
	  else {
	    bool s((*b)(0)==(*c)(0));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*b)(0))) {
		m_ci.push_back(CInfo((*b)(0),i,1,0));
		m_cjk.push_back(CInfo(i,(*c)(1),0,0));
	      }
	    if (!s) m_ci.push_back(CInfo((*c)(0),(*a)(1),0,0));
	  }
	  if ((*b)(0)!=(*c)(1)) {
	    m_ci.push_back(CInfo((*b)(0),(*c)(1),1,1));
	    m_cjk.push_back(CInfo((*c)(0),(*a)(1),1,1));
	  }
	  else {
	    bool s((*a)(1)==(*c)(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*a)(1))) {
		m_ci.push_back(CInfo(i,(*a)(1),0,1));
		m_cjk.push_back(CInfo((*c)(0),i,1,1));
	      }
	    if (!s) m_ci.push_back(CInfo((*b)(0),(*c)(1),1,1));
	  }
	  break;
	}
	default: THROW(fatal_error,"Invalid call");
	}
	break;
      }
      default: THROW(fatal_error,"Invalid call");
      }
      m_stat=m_cjk.size();
      return m_stat;
    }

  };// end of class SF_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(SF_C_Getter,"S-F",Color_Calculator,Vertex_Key);

Color_Calculator *SF_C_Getter::operator()(const Vertex_Key &key) const
{
  return new SF_Calculator(key);
}

void SF_C_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"adjoint (subtraction)";
}

