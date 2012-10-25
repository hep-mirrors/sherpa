#ifndef PHASIC_Selectors_Decay2_Selector_H
#define PHASIC_Selectors_Decay2_Selector_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {

  class Decay2_Selector: public Selector_Base,
			 public ATOOLS::Tag_Replacer {
  private:

    std::vector<std::vector<int> > m_ids[2];

    ATOOLS::Vec4D_Vector m_p[2];

    double m_min, m_max;

    ATOOLS::Algebra_Interpreter m_calc;

  public:

    Decay2_Selector(const Selector_Key &key);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);

    void BuildCuts(Cut_Data *) {}

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

  };

}

#endif

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace PHASIC;
using namespace ATOOLS;

Decay2_Selector::Decay2_Selector(const Selector_Key &key):
  Selector_Base("Decay2_Selector")
{
  if (key.m_key.length()<7 || key[0].size()<3)
    THROW(fatal_error,"Invalid syntax");
  std::string tag(key.m_key.substr(7));
  tag.erase(tag.length()-1,1);
  DEBUG_FUNC(tag);
  long int kf1(ToType<long int>(key[0][0]));
  long int kf2(ToType<long int>(key[0][1]));
  Flavour fl1(abs(kf1),kf1<0), fl2(abs(kf2),kf2<0);
  DecayInfo_Vector decs
    (key.p_proc->Process()->Info().m_fi.GetDecayInfos());
  for (size_t i(0);i<decs.size();++i) {
    if (decs[i]->m_fl==fl1) {
      m_ids[0].push_back(ID(decs[i]->m_id));
      if (m_ids[0].size()>1 &&
	  m_ids[0].front().size()!=m_ids[0].back().size())
	THROW(fatal_error,"Varying multiplicity");
      msg_Debugging()<<"adding "<<m_ids[0].back()<<"\n";
    }
    if (decs[i]->m_fl==fl2) {
      m_ids[1].push_back(ID(decs[i]->m_id));
      if (m_ids[1].size()>1 &&
	  m_ids[1].front().size()!=m_ids[1].back().size())
	THROW(fatal_error,"Varying multiplicity");
      msg_Debugging()<<"adding "<<m_ids[1].back()<<"\n";
    }
  }
  if (m_ids[0].empty() || m_ids[1].empty())
    THROW(fatal_error,"No such flavour");
  m_p[0].resize(m_ids[0].back().size());
  m_p[1].resize(m_ids[1].back().size());
  for (size_t i(0);i<m_p[0].size();++i) 
    m_calc.AddTag("p1["+ToString(i)+"]",ToString(Vec4D()));
  for (size_t i(0);i<m_p[1].size();++i) 
    m_calc.AddTag("p2["+ToString(i)+"]",ToString(Vec4D()));
  m_calc.SetTagReplacer(this);
  m_calc.Interprete(tag);
  if (msg_LevelIsDebugging()) m_calc.PrintEquation();
  m_min=ToType<double>(key[0][2]);
  m_max=ToType<double>(key[0][3]);
  msg_Debugging()<<"m_min = "<<m_min
		 <<", m_max = "<<m_max<<"\n";
  m_sel_log = new Selector_Log(m_name);
}

bool Decay2_Selector::Trigger(const Vec4D_Vector &p)
{
  DEBUG_FUNC("");
  for (size_t j(0);j<m_ids[0].size();++j) {
    for (size_t i(0);i<m_ids[0][j].size();++i) m_p[0][i]=p[m_ids[0][j][i]];
    for (size_t l(0);l<m_ids[1].size();++l) {
      for (size_t k(0);k<m_ids[1][l].size();++k) m_p[1][k]=p[m_ids[1][l][k]];
      double value(m_calc.Calculate()->Get<double>());
      msg_Debugging()<<m_ids[0][j]<<","<<m_ids[1][l]<<" -> "<<value<<"\n";
      if (value<m_min || value>m_max) return !m_sel_log->Hit(1);
    }
  }
  return !m_sel_log->Hit(0);
}

std::string Decay2_Selector::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Decay2_Selector::ReplaceTags(Term *term) const
{
  if (term->Id()>=200) term->Set(m_p[1][term->Id()-200]);
  else if (term->Id()>=100) term->Set(m_p[0][term->Id()-100]);
  return term;
}

void Decay2_Selector::AssignId(Term *term)
{
  if (term->Tag().find("p1")==0) {
    term->SetId(100+ToType<int>
		(term->Tag().substr
		 (3,term->Tag().length()-4)));
  }
  else if (term->Tag().find("p2")==0) {
    term->SetId(200+ToType<int>
		(term->Tag().substr
		 (3,term->Tag().length()-4)));
  }
}

namespace PHASIC {

DECLARE_ND_GETTER(Decay2_Selector_Getter,"Decay2",
		  Selector_Base,Selector_Key,true);

Selector_Base *Decay2_Selector_Getter::operator()
  (const Selector_Key &key) const
{
  Decay2_Selector *msel(new Decay2_Selector(key));
  return msel;
}

void Decay2_Selector_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"Decay2 kf1 kf2 min max"; 
}

}
