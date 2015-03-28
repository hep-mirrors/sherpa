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
      m_ti=key.Fl(0).StrongCharge();
      m_tj=key.Fl(1).StrongCharge();
      m_tk=key.p_k->Flav().StrongCharge();
    }

    std::string Label() const
    {
      return "S-T";
    }

    bool Evaluate(const CObject_Vector &j)
    {
      if (j.size()==2) {
      m_ci.clear();
      m_cjk.clear();
      m_stat=true;
      if (m_ti==3 || m_tj==3) {
	switch (m_tk) {
	case -3:
	  if ((*j[0])(0)!=(*j[1])(1)) {
	    m_ci.push_back(CInfo((*j[0])(0),0,0,0,-3.0));
	    m_cjk.push_back(CInfo(0,(*j[1])(1),1,0));
	  }
	  else {
	    for (size_t i(s_cimin);i<=s_cimax;++i) {
	      m_ci.push_back(CInfo(i,0,0,0));
	      m_cjk.push_back(CInfo(0,i,1,0));
	    }
	    m_ci.push_back(CInfo((*j[0])(0),0,0,0,-3.0));
	  }
	  break;
	case 3:
	  m_cjk.push_back(CInfo((*j[0])(0),0,0,0));
	  m_cjk.push_back(CInfo((*j[1])(0),0,0,1)); 
	  m_ci.push_back(CInfo((*j[1])(0),0,0,0));
	  m_ci.push_back(CInfo((*j[0])(0),0,0,1,-3.0));
	  break;
	case 8:
	  if ((*j[0])(0)!=(*j[1])(1)) {
	    m_ci.push_back(CInfo((*j[1])(0),0,0,0));
	    m_cjk.push_back(CInfo((*j[0])(0),(*j[1])(1),0,0));
	  }
	  else {
	    bool s((*j[1])(0)==(*j[1])(1));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*j[1])(0))) {
		m_ci.push_back(CInfo(i,0,0,0));
		m_cjk.push_back(CInfo((*j[1])(0),i,1,0));
	      }
	    if (!s) {
	      m_ci.push_back(CInfo((*j[1])(0),0,0,1));
	      m_cjk.push_back(CInfo((*j[0])(0),(*j[1])(1),0,1));
	    }
	  }
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
      }
      else if (m_ti==-3 || m_tj==-3) {
	switch (m_tk) {
	case 3:
	  if ((*j[0])(1)!=(*j[1])(0)) {
	    m_ci.push_back(CInfo(0,(*j[0])(1),1,0,-3.0));
	    m_cjk.push_back(CInfo((*j[1])(0),0,0,0));
	  }
	  else {
	    for (size_t i(s_cimin);i<=s_cimax;++i) {
	      m_ci.push_back(CInfo(0,i,1,0));
	      m_cjk.push_back(CInfo(i,0,0,0));
	    }
	    m_ci.push_back(CInfo(0,(*j[0])(1),1,0,-3.0));
	  }
	  break;
	case -3:
	  m_cjk.push_back(CInfo(0,(*j[0])(1),1,0));
	  m_cjk.push_back(CInfo(0,(*j[1])(1),1,1));
	  m_ci.push_back(CInfo(0,(*j[1])(1),1,0));
	  m_ci.push_back(CInfo(0,(*j[0])(1),1,1,-3.0));
	  break;
	case 8:
	  if ((*j[0])(1)!=(*j[1])(0)) {
	    m_ci.push_back(CInfo(0,(*j[1])(1),1,0));
	    m_cjk.push_back(CInfo((*j[1])(0),(*j[0])(1),1,0));
	  }
	  else {
	    bool s((*j[1])(1)==(*j[1])(0));
	    for (size_t i(s_cimin);i<=s_cimax;++i)
	      if (!(s && (int)i==(*j[1])(1))) {
		m_ci.push_back(CInfo(0,i,1,0));
		m_cjk.push_back(CInfo(i,(*j[1])(1),0,0));
	      }
	    if (!s) {
	      m_ci.push_back(CInfo(0,(*j[1])(1),1,1));
	      m_cjk.push_back(CInfo((*j[1])(0),(*j[0])(1),1,1));
	    }
	  }
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
      }
      else {
	THROW(fatal_error,"Invalid call");
      }
      return m_stat;
      }
      m_ci.clear();
      m_cjk.clear();
      m_stat=p_cc->Evaluate(j);
      if (!m_stat) return false;
      switch (m_ti) {
      case -3:
      case 3:
	m_s=(*j[1])(0)==(*j[1])(1);
	m_mj=(*j[1])(0)==(*j[2])(1);
	m_mi=(*j[1])(1)==(*j[2])(0);
	switch (m_tk) {
	case -3:
	  m_stat=m_mj||m_s;
	  if (m_mj) m_cjk.push_back(CInfo(0,(*j[1])(1),1,0));
	  if (m_s) m_cjk.push_back(CInfo(0,(*j[2])(1),1,0,-3.0));
	  break;
	case 3:
	  m_stat=m_mi||m_s;
	  if (m_mi) m_cjk.push_back(CInfo((*j[1])(0),0,0,0));
	  if (m_s) m_cjk.push_back(CInfo((*j[2])(0),0,0,0,-3.0));
	  break;
	case 8:
	  m_stat=(m_mi||m_mj)&&!(m_mi&&m_mj&&m_s);
	  if (m_mj) m_cjk.push_back(CInfo((*j[2])(0),(*j[1])(1),1,0));
	  if (m_mi) m_cjk.push_back(CInfo((*j[1])(0),(*j[2])(1),0,0));
	  break;
	default: THROW(fatal_error,"Invalid call");
	}
	if (m_stat) {
	  if (m_ti>0) m_ci.push_back(CInfo((*j[0])(0),0,0,0));
	  if (m_ti<0) m_ci.push_back(CInfo(0,(*j[0])(1),1,0));
	}
	break;
      case 8:
	m_s=(*j[0])(0)==(*j[0])(1);
	  m_mj=(*j[0])(1)==(*j[2])(0);
	  m_mi=(*j[0])(0)==(*j[2])(1);
	  switch (m_tk) {
	  case 3:
	    m_stat=m_mj||m_s;
	    if (m_mj) m_cjk.push_back(CInfo((*j[0])(0),0,0,0));
	    if (m_s) m_cjk.push_back(CInfo((*j[2])(0),0,0,0,-3.0));
	    break;
	  case -3:
	    m_stat=m_mi||m_s;
	    if (m_mi) m_cjk.push_back(CInfo(0,(*j[0])(1),1,0));
	    if (m_s) m_cjk.push_back(CInfo(0,(*j[2])(1),1,0,-3.0));
	    break;
	  case 8:
	    m_stat=(m_mi||m_mj)&&!(m_mi&&m_mj&&m_s);
	    if (m_mj) m_cjk.push_back(CInfo((*j[0])(0),(*j[2])(1),0,0));
	    if (m_mi) m_cjk.push_back(CInfo((*j[2])(0),(*j[0])(1),1,0));
	    break;
	  default: THROW(fatal_error,"Invalid call");
	  }
	  if (m_stat) {
	    if (m_tj>0) m_ci.push_back(CInfo((*j[1])(0),0,0,0));
	    if (m_tj<0) m_ci.push_back(CInfo(0,(*j[1])(1),1,0));
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

DECLARE_GETTER(ST_Calculator,"S-T",
	       Color_Calculator,Vertex_Key);

Color_Calculator *ATOOLS::Getter
<Color_Calculator,Vertex_Key,ST_Calculator>::
operator()(const Vertex_Key &key) const
{
  return new ST_Calculator(key);
}

void ATOOLS::Getter<Color_Calculator,Vertex_Key,ST_Calculator>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"fundamental (subtraction)";
}
