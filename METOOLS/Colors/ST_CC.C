#include "METOOLS/Explicit/Dipole_Color.H"

#include "METOOLS/Explicit/Vertex.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

namespace METOOLS {

  class ST_Calculator: public Dipole_Color {
  private:

    int m_ti, m_tj, m_tk, m_mi, m_mj, m_s;

  public:

    inline ST_Calculator(const Vertex_Key &key): 
      Dipole_Color(key)
    {
      m_cpl=p_cc->Coupling();
      if (key.p_c->Flav().StrongCharge()==8)
	THROW(fatal_error,"Invalid call");
      m_ti=key.FlA().StrongCharge();
      m_tj=key.FlB().StrongCharge();
      m_tk=key.p_k->Flav().StrongCharge();
    }

    std::string Label() const
    {
      return "S-T";
    }

    bool Evaluate(const CObject *a,const CObject *b)
    {
      m_ci.clear();
      m_cjk.clear();
      m_stat=true;
      switch (m_ti) {
      case 3:
	switch (m_tk) {
	case -3:
	  if ((*a)(0)!=(*b)(1)) {
	    m_ci.push_back(CInfo((*a)(0),0,0,0,-3.0));
	    m_cjk.push_back(CInfo(0,(*b)(1),1,0));
	  }
	  else {
	    for (size_t i(s_cimin);i<=s_cimax;++i) {
	      m_ci.push_back(CInfo(i,0,0,0));
	      m_cjk.push_back(CInfo(0,i,1,0));
	    }
	    m_ci.push_back(CInfo((*a)(0),0,0,0,-3.0));
	  }
	  break;
	case 3:
	  m_cjk.push_back(CInfo((*a)(0),0,0,0));
	  m_cjk.push_back(CInfo((*b)(0),0,0,1)); 
	  m_ci.push_back(CInfo((*b)(0),0,0,0));
	  m_ci.push_back(CInfo((*a)(0),0,0,1,-3.0));
	  break;
	case 8:
	  if ((*a)(0)!=(*b)(1)) {
	    m_ci.push_back(CInfo((*b)(0),0,0,0));
	    m_cjk.push_back(CInfo((*a)(0),(*b)(1),0,0));
	  }
	  else {
	    bool s((*b)(0)==(*b)(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*b)(0))) {
		m_ci.push_back(CInfo(i,0,0,0));
		m_cjk.push_back(CInfo((*b)(0),i,1,0));
	      }
	    if (!s) {
	      m_ci.push_back(CInfo((*b)(0),0,0,1));
	      m_cjk.push_back(CInfo((*a)(0),(*b)(1),0,1));
	    }
	  }
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	break;
      case 8:
	switch (m_tj) {
	case -3:
	  switch (m_tk) {
	  case 3:
	    if ((*a)(1)!=(*b)(0)) {
	      m_ci.push_back(CInfo(0,(*a)(1),1,0,-3.0));
	      m_cjk.push_back(CInfo((*b)(0),0,0,0));
	    }
	    else {
	      for (size_t i(s_cimin);i<=s_cimax;++i) {
		m_ci.push_back(CInfo(0,i,1,0));
		m_cjk.push_back(CInfo(i,0,0,0));
	      }
	      m_ci.push_back(CInfo(0,(*a)(1),1,0,-3.0));
	    }
	    break;
	  case -3:
	    m_cjk.push_back(CInfo(0,(*a)(1),1,0));
	    m_cjk.push_back(CInfo(0,(*b)(1),1,1));
	    m_ci.push_back(CInfo(0,(*b)(1),1,0));
	    m_ci.push_back(CInfo(0,(*a)(1),1,1,-3.0));
	    break;
	  case 8:
	    if ((*a)(1)!=(*b)(0)) {
	      m_ci.push_back(CInfo(0,(*b)(1),1,0));
	      m_cjk.push_back(CInfo((*b)(0),(*a)(1),1,0));
	    }
	    else {
	      bool s((*b)(1)==(*b)(0));
	      for (size_t i(s_cimin);i<=s_cimax;++i)
		if (!(s && (int)i==(*b)(1))) {
		  m_ci.push_back(CInfo(0,i,1,0));
		  m_cjk.push_back(CInfo(i,(*b)(1),0,0));
		}
	      if (!s) {
		m_ci.push_back(CInfo(0,(*b)(1),1,1));
		m_cjk.push_back(CInfo((*b)(0),(*a)(1),1,1));
	      }
	    }
	    break;
	  default: THROW(fatal_error,"Invalid call");
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
      case 3:
	m_s=(*b)(0)==(*b)(1);
	m_mj=(*b)(0)==(*c)(1);
	m_mi=(*b)(1)==(*c)(0);
	switch (m_tk) {
	case -3:
	  m_stat=m_mj||m_s;
	  if (m_mj) m_cjk.push_back(CInfo(0,(*b)(1),1,0));
	  if (m_s) m_cjk.push_back(CInfo(0,(*c)(1),1,0,-3.0));
	  break;
	case 3:
	  m_stat=m_mi||m_s;
	  if (m_mi) m_cjk.push_back(CInfo((*b)(0),0,0,0));
	  if (m_s) m_cjk.push_back(CInfo((*c)(0),0,0,0,-3.0));
	  break;
	case 8:
	  m_stat=(m_mi||m_mj)&&!(m_mi&&m_mj&&m_s);
	  if (m_mj) m_cjk.push_back(CInfo((*c)(0),(*b)(1),1,0));
	  if (m_mi) m_cjk.push_back(CInfo((*b)(0),(*c)(1),0,0));
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	if (m_stat) m_ci.push_back(CInfo((*a)(0),0,0,0));
	break;
      case 8:
	m_s=(*a)(0)==(*a)(1);
	switch (m_tj) {
	case -3:
	  m_mj=(*a)(1)==(*c)(0);
	  m_mi=(*a)(0)==(*c)(1);
	  switch (m_tk) {
	  case 3:
	    m_stat=m_mj||m_s;
	    if (m_mj) m_cjk.push_back(CInfo((*a)(0),0,0,0));
	    if (m_s) m_cjk.push_back(CInfo((*c)(0),0,0,0,-3.0));
	    break;
	  case -3:
	    m_stat=m_mi||m_s;
	    if (m_mi) m_cjk.push_back(CInfo(0,(*a)(1),1,0));
	    if (m_s) m_cjk.push_back(CInfo(0,(*c)(1),1,0,-3.0));
	    break;
	  case 8:
	    m_stat=(m_mi||m_mj)&&!(m_mi&&m_mj&&m_s);
	    if (m_mj) m_cjk.push_back(CInfo((*a)(0),(*c)(1),0,0));
	    if (m_mi) m_cjk.push_back(CInfo((*c)(0),(*a)(1),1,0));
	    break;
	  default: THROW(fatal_error,"Invalid call");
	  }
	  if (m_stat) m_ci.push_back(CInfo(0,(*b)(1),1,0));
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	break;
      default: THROW(fatal_error,"Invalid call");
      }
      return m_stat;
    }

  };// end of class ST_Calculator

}// end of namespace METOOLS

using namespace METOOLS;
using namespace ATOOLS;

DECLARE_GETTER(ST_C_Getter,"S-T",Color_Calculator,Vertex_Key);

Color_Calculator *ST_C_Getter::operator()(const Vertex_Key &key) const
{
  return new ST_Calculator(key);
}

void ST_C_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"fundamental (subtraction)";
}

