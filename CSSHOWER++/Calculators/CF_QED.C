#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"

namespace CSSHOWER {
  
  class CF_QED: public SF_Coupling {
  protected:

    ATOOLS::Function_Base *p_cpl;
    ATOOLS::Flavour m_cfl;

    double m_q;

  public:

    inline CF_QED(const SF_Key &key):
      SF_Coupling(key), m_cfl(key.p_v->in[0].Bar())
    {
      if (key.m_type==cstp::IF || key.m_type==cstp::II)
	m_cfl=key.p_v->in[key.m_mode==0?1:2];
      m_q=ATOOLS::dabs(m_cfl.IntCharge()?m_cfl.Charge():
		       key.p_v->in[key.m_mode==0?2:1].Charge());
      if (m_q==0.0) THROW(fatal_error,"Internal error");
    }

    bool SetCoupling(MODEL::Model_Base *md,
		     const double &k0sqi,const double &k0sqf,
		     const double &isfac,const double &fsfac);
    double Coupling(const double &scale,const int pol);
    bool AllowSpec(const ATOOLS::Flavour &fl);

  };

}

using namespace CSSHOWER;
using namespace MODEL;
using namespace ATOOLS;

bool CF_QED::SetCoupling(MODEL::Model_Base *md,
			 const double &k0sqi,const double &k0sqf,
			 const double &isfac,const double &fsfac)
{
  p_cpl=md->GetScalarFunction("alpha_QED");
  m_cplfac=1.0;
  m_cplmax.push_back((*p_cpl)(sqr(rpa->gen.Ecms()))*m_q);
  m_cplmax.push_back(0.0);
  return true;
}

double CF_QED::Coupling(const double &scale,const int pol)
{
  if (pol!=0) return 0.0;
  if (scale<0.0) return
    m_cplmax.front()*m_q*dabs(p_lf->FlSpec().Charge());
  double scl(CplFac(scale)*scale);
  return (*p_cpl)(scl)*m_q*dabs(p_lf->FlSpec().Charge());
}

bool CF_QED::AllowSpec(const ATOOLS::Flavour &fl) 
{
  if (!fl.Strong() && fl.Mass()>10.0) return false;
  if (m_cfl.IntCharge()==0) return fl.Charge();
  
  switch (m_type) {
  case cstp::FF:
  case cstp::II:
    return fl.IntCharge()*m_cfl.IntCharge()<0;
  case cstp::FI:
  case cstp::IF:
    return fl.IntCharge()*m_cfl.IntCharge()>0;
  default:
    return false;
  }
}

namespace CSSHOWER {

DECLARE_CPL_GETTER(CF_QED_Getter);

SF_Coupling *CF_QED_Getter::operator()
  (const Parameter_Type &args) const
{
  return new CF_QED(args);
}

void CF_QED_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{
  str<<"electromagnetic coupling";
}

}

DECLARE_GETTER(CF_QED_Getter,"SF_QED_Fill",
	       void,SFC_Filler_Key);

void *ATOOLS::Getter<void,SFC_Filler_Key,CF_QED_Getter>::
operator()(const SFC_Filler_Key &key) const
{
  DEBUG_FUNC("model = "<<key.p_md->Name());
  const Vertex_Table *vtab(key.p_md->VertexTable());
  for (Vertex_Table::const_iterator
	 vlit=vtab->begin();vlit!=vtab->end();++vlit) {
    for (Vertex_List::const_iterator 
	   vit=vlit->second.begin();vit!=vlit->second.end();++vit) {
      Single_Vertex *v(*vit);
      if (v->NLegs()>3) continue;
      if (!((v->in[0].IsPhoton() && (v->in[1].IsFermion()||v->in[1].IsScalar()) && v->in[1].Charge()) ||
	    (v->in[1].IsPhoton() && (v->in[2].IsFermion()||v->in[2].IsScalar()) && v->in[2].Charge()) ||
	    (v->in[2].IsPhoton() && (v->in[0].IsFermion()||v->in[0].IsScalar()) && v->in[0].Charge()))) continue;
      msg_Debugging()<<"Add "<<v->in[0].Bar()<<" -> "<<v->in[1]<<" "<<v->in[2]<<" {\n";
      std::string atag("{"+v->in[0].Bar().IDName()+"}");
      std::string btag("{"+v->in[1].IDName()+"}");
      std::string ctag("{"+v->in[2].IDName()+"}");
      key.p_gets->push_back(new CF_QED_Getter(atag+btag+ctag));
    }
  }
  return NULL;
}

void ATOOLS::Getter<void,SFC_Filler_Key,CF_QED_Getter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"qed coupling filler";
}


