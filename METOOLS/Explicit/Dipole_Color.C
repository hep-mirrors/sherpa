#include "METOOLS/Explicit/Dipole_Color.H"

#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Explicit/Vertex.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"

#define COMPILE__Getter_Function
#define PARAMETER_TYPE METOOLS::Vertex_Key
#define OBJECT_TYPE METOOLS::Dipole_Color
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"

using namespace METOOLS;

Dipole_Color::Dipole_Color(const Vertex_Key &key): 
  Color_Calculator(key), p_dinfo(key.p_dinfo), p_kt(key.p_kt)
{
  std::string ctag(ToString(key.p_mv->Color[key.m_n].PID()));
  key.m_d=1;
  p_cc = CC_Getter::GetObject(ctag,key);
  key.m_d=0;
  if (p_cc==NULL) {
    msg_Info()<<*key.p_mv<<std::endl;
    THROW(fatal_error,"Color calculator not implemented '"+
	  ctag+"'");
  }
}

Dipole_Color::~Dipole_Color()
{
  delete p_cc;
}

void Dipole_Color::AddJ(CObject *const j)
{
  p_cc->AddJ(j);
}

void Dipole_Color::AddJI(CObject *const j,const int t) const
{
  for (std::vector<CInfo>::const_iterator
	 c(m_ci.begin());c<m_ci.end();++c)
    if (c->m_t==t) {
      CObject *cc(j->Copy());
      cc->SetS(1+c->m_t);
      (*cc)(0)=c->m_cr;
      (*cc)(1)=c->m_ca;
      if (c->m_s!=1.0) cc->Divide(c->m_s);
      if (c->m_i) cc->Invert();
      p_v->AddJ(cc);
    }
}

void Dipole_Color::AddJJK(CObject *const j) const
{
  for (std::vector<CInfo>::const_iterator
	 c(m_cjk.begin());c<m_cjk.end();++c) {
    CObject *cc(j->Copy());
    cc->SetS(1+c->m_t);
    (*cc)(0)=c->m_cr;
    (*cc)(1)=c->m_ca;
    if (c->m_s!=1.0) cc->Divide(c->m_s);
    if (c->m_i) cc->Invert();
    p_kt->AddJ(cc);
  }
}

namespace METOOLS {

  std::ostream &operator<<(std::ostream &str,
			   const Dipole_Color::CInfo &c)
  {
    return str<<'{'<<c.m_cr<<','<<c.m_ca<<'|'
	      <<c.m_i<<','<<c.m_s<<';'<<c.m_t<<'}';
  }

}
