#ifndef PHASIC_Selectors_Decay_Selector_H
#define PHASIC_Selectors_Decay_Selector_H

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {

  class Decay_Selector: public Selector_Base,
			public ATOOLS::Tag_Replacer {
  private:

    std::vector<std::vector<int> > m_ids;

    ATOOLS::Vec4D_Vector m_p;

    double m_min, m_max;

    ATOOLS::Algebra_Interpreter m_calc;

  public:

    Decay_Selector(const Selector_Key &key);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);
    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
    { return Trigger(p); }
    bool JetTrigger(const ATOOLS::Vec4D_Vector &p,
		    ATOOLS::NLO_subevtlist *const subs)
    { return true; }

    void BuildCuts(Cut_Data *) {}

    std::string   ReplaceTags(std::string &expr) const;    
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;    
    
    void AssignId(ATOOLS::Term *term);

  };

  class DecayMass_Selector: public Selector_Base {
  private:

    std::vector<std::vector<int> > m_ids;

    double m_min, m_max;

  public:

    DecayMass_Selector(const Selector_Key &key);

    bool Trigger(const ATOOLS::Vec4D_Vector &p);
    bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
    { return Trigger(p); }
    bool JetTrigger(const ATOOLS::Vec4D_Vector &p,
		    ATOOLS::NLO_subevtlist *const subs)
    { return true; }

    void BuildCuts(Cut_Data *);

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

Decay_Selector::Decay_Selector(const Selector_Key &key):
  Selector_Base("Decay_Selector")
{
  if (key.m_key.length()<7 || key[0].size()<3)
    THROW(fatal_error,"Invalid syntax");
  std::string tag(key.m_key.substr(6));
  tag.erase(tag.length()-1,1);
  DEBUG_FUNC(tag);
  long int kf(ToType<long int>(key[0][0]));
  Flavour fl(abs(kf),kf<0);
  DecayInfo_Vector decs
    (key.p_proc->Process()->Info().m_fi.GetDecayInfos());
  for (size_t i(0);i<decs.size();++i)
    if (decs[i]->m_fl==fl) {
      m_ids.push_back(ID(decs[i]->m_id));
      if (m_ids.size()>1 &&
	  m_ids.front().size()!=m_ids.back().size())
	THROW(fatal_error,"Varying multiplicity");
      msg_Debugging()<<"adding "<<m_ids.back()<<"\n";
    }
  if (m_ids.empty()) THROW(fatal_error,"No such flavour");
  m_p.resize(m_ids.back().size());
  for (size_t i(0);i<m_p.size();++i) 
    m_calc.AddTag("p["+ToString(i)+"]",ToString(Vec4D()));
  m_calc.SetTagReplacer(this);
  m_calc.Interprete(tag);
  if (msg_LevelIsDebugging()) m_calc.PrintEquation();
  m_min=ToType<double>(key[0][1]);
  m_max=ToType<double>(key[0][2]);
  msg_Debugging()<<"m_min = "<<m_min
		 <<", m_max = "<<m_max<<"\n";
  m_sel_log = new Selector_Log(m_name);
}

bool Decay_Selector::Trigger(const Vec4D_Vector &p)
{
  DEBUG_FUNC("");
  for (size_t j(0);j<m_ids.size();++j) {
    for (size_t i(0);i<m_ids[j].size();++i) m_p[i]=p[m_ids[j][i]];
    double value(m_calc.Calculate()->Get<double>());
    msg_Debugging()<<m_ids[j]<<" -> "<<value<<"\n";
    if (value<m_min || value>m_max) return !m_sel_log->Hit(1);
  }
  return !m_sel_log->Hit(0);
}

std::string Decay_Selector::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Decay_Selector::ReplaceTags(Term *term) const
{
  term->Set(m_p[term->Id()]);
  return term;
}

void Decay_Selector::AssignId(Term *term)
{
  term->SetId(ToType<int>
	      (term->Tag().substr
	       (2,term->Tag().length()-3)));
}

namespace PHASIC {

DECLARE_ND_GETTER(Decay_Selector_Getter,"Decay",
		  Selector_Base,Selector_Key,true);

Selector_Base *Decay_Selector_Getter::operator()
  (const Selector_Key &key) const
{
  Decay_Selector *msel(new Decay_Selector(key));
  return msel;
}

void Decay_Selector_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"Decay kf min max"; 
}

}

DecayMass_Selector::DecayMass_Selector(const Selector_Key &key):
  Selector_Base("DecayMass_Selector")
{
  DEBUG_FUNC(key.m_key);
  long int kf(ToType<long int>(key[0][0]));
  Flavour fl(abs(kf),kf<0);
  DecayInfo_Vector decs
    (key.p_proc->Process()->Info().m_fi.GetDecayInfos());
  for (size_t i(0);i<decs.size();++i)
    if (decs[i]->m_fl==fl) {
      m_ids.push_back(ID(decs[i]->m_id));
      if (m_ids.size()>1 &&
	  m_ids.front().size()!=m_ids.back().size())
	THROW(fatal_error,"Varying multiplicity");
      msg_Debugging()<<"adding "<<m_ids.back()<<"\n";
    }
  if (m_ids.empty()) THROW(fatal_error,"No such flavour");
  m_min=ToType<double>(key[0][1]);
  m_max=ToType<double>(key[0][2]);
  msg_Debugging()<<"m_min = "<<m_min
		 <<", m_max = "<<m_max<<"\n";
  m_sel_log = new Selector_Log(m_name);
}

bool DecayMass_Selector::Trigger(const Vec4D_Vector &p)
{
  DEBUG_FUNC("");
  Vec4D sum;
  for (size_t j(0);j<m_ids.size();++j) {
    for (size_t i(0);i<m_ids[j].size();++i) sum+=p[m_ids[j][i]];
    double value(sum.Mass());
    msg_Debugging()<<m_ids[j]<<" -> "<<value<<"\n";
    if (value<m_min || value>m_max) return !m_sel_log->Hit(1);
  }
  return !m_sel_log->Hit(0);
}

void DecayMass_Selector::BuildCuts(Cut_Data *cuts)
{
  for (size_t j(0);j<m_ids.size();++j) {
    if (m_ids[j].size()==2) {
      cuts->scut[m_ids[j][0]][m_ids[j][1]]=
	cuts->scut[m_ids[j][1]][m_ids[j][0]]=
	Max(cuts->scut[m_ids[j][0]][m_ids[j][1]],sqr(m_min));
    }
    std::string id;
    for (size_t i(0);i<m_ids[j].size();++i) id+=ToString(m_ids[j][i]);
    double scut(cuts->Getscut(id));
    cuts->Setscut(id,Max(scut,sqr(m_min)));
  }
}

namespace PHASIC {

DECLARE_ND_GETTER(DecayMass_Selector_Getter,"DecayMass",
		  Selector_Base,Selector_Key,true);

Selector_Base *DecayMass_Selector_Getter::operator()
  (const Selector_Key &key) const
{
  DecayMass_Selector *msel(new DecayMass_Selector(key));
  return msel;
}

void DecayMass_Selector_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"DecayMass kf min max"; 
}

}
