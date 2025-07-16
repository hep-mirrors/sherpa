#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"

#include <cassert>

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

    bool Trigger(ATOOLS::Selector_List &p);

    void BuildCuts(Cut_Data *) {}

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  };

  class Decay2_Selector: public Selector_Base,
                         public ATOOLS::Tag_Replacer {
  private:

    std::vector<std::vector<int> > m_ids[2];

    ATOOLS::Vec4D_Vector m_p[2];

    double m_min, m_max;

    ATOOLS::Algebra_Interpreter m_calc;

  public:

    Decay2_Selector(const Selector_Key &key);

    bool Trigger(ATOOLS::Selector_List &p);

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

    bool Trigger(ATOOLS::Selector_List &p);

    void BuildCuts(Cut_Data *);

  };

}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

Decay_Selector::Decay_Selector(const Selector_Key &key):
  Selector_Base("Decay_Selector",key.p_proc)
{
  Scoped_Settings s{ key.m_settings };
  s.SetInterpreterEnabled(false);
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters[0].size() < 7 || parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  std::string tag(parameters[0].substr(6));
  tag.pop_back();
  DEBUG_FUNC(tag);
  const auto kf = s.Interprete<long int>(parameters[1]);
  Flavour fl = Flavour((kf_code)std::abs(kf),kf<0);
  DecayInfo_Vector decs(p_proc->Info().m_fi.GetDecayInfos());
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
  m_min = s.Interprete<double>(parameters[2]);
  m_max = s.Interprete<double>(parameters[3]);
  msg_Debugging()<<"m_min = "<<m_min
		 <<", m_max = "<<m_max<<"\n";
}

bool Decay_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC("");
  for (size_t j(0);j<m_ids.size();++j) {
    for (size_t i(0);i<m_ids[j].size();++i) m_p[i]=sl[m_ids[j][i]].Momentum();
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

DECLARE_GETTER(Decay_Selector,"Decay",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Decay_Selector>::
operator()(const Selector_Key &key) const
{
  Decay_Selector *msel(new Decay_Selector(key));
  return msel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Decay_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[Decay, kf, min, max]";
}


Decay2_Selector::Decay2_Selector(const Selector_Key &key):
  Selector_Base("Decay2_Selector",key.p_proc)
{
  Scoped_Settings s{ key.m_settings };
  s.SetInterpreterEnabled(false);
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters[0].size() < 7 || parameters.size() != 5)
    THROW(critical_error, "Invalid syntax");
  std::string tag(parameters[0].substr(7));
  tag.pop_back();
  DEBUG_FUNC(tag);
  const auto kf1 = s.Interprete<long int>(parameters[1]);
  const auto kf2 = s.Interprete<long int>(parameters[2]);
  Flavour fl1 = Flavour((kf_code)std::abs(kf1),kf1<0);
  Flavour fl2 = Flavour((kf_code)std::abs(kf2),kf2<0);
  DecayInfo_Vector decs(p_proc->Info().m_fi.GetDecayInfos());
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
  m_min = s.Interprete<double>(parameters[3]);
  m_max = s.Interprete<double>(parameters[4]);
  msg_Debugging()<<"m_min = "<<m_min
                 <<", m_max = "<<m_max<<"\n";
}

bool Decay2_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC("");
  for (size_t j(0);j<m_ids[0].size();++j) {
    for (size_t i(0);i<m_ids[0][j].size();++i) m_p[0][i]=sl[m_ids[0][j][i]].Momentum();
    for (size_t l(0);l<m_ids[1].size();++l) {
      for (size_t k(0);k<m_ids[1][l].size();++k) m_p[1][k]=sl[m_ids[1][l][k]].Momentum();
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

DECLARE_GETTER(Decay2_Selector,"Decay2",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Decay2_Selector>::
operator()(const Selector_Key &key) const
{
  Decay2_Selector *msel(new Decay2_Selector(key));
  return msel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Decay2_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[Decay2, kf1, kf2, min, max]";
}


DecayMass_Selector::DecayMass_Selector(const Selector_Key &key):
  Selector_Base("DecayMass_Selector",key.p_proc)
{
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DecayMass");
  if (parameters.size() != 4)
    THROW(critical_error, "Invalid syntax");
  const auto kf = s.Interprete<int>(parameters[1]);
  Flavour fl = Flavour((kf_code)abs(kf),kf<0);
  DecayInfo_Vector decs(p_proc->Info().m_fi.GetDecayInfos());
  for (size_t i(0);i<decs.size();++i)
    if (decs[i]->m_fl==fl) {
      m_ids.push_back(ID(decs[i]->m_id));
      if (m_ids.size()>1 &&
	  m_ids.front().size()!=m_ids.back().size())
	THROW(fatal_error,"Varying multiplicity");
      msg_Debugging()<<"adding "<<m_ids.back()<<"\n";
    }
  if (m_ids.empty()) THROW(fatal_error,"No such flavour");
  m_min = s.Interprete<double>(parameters[2]);
  m_max = s.Interprete<double>(parameters[3]);
  msg_Debugging()<<"m_min = "<<m_min
		 <<", m_max = "<<m_max<<"\n";
}

bool DecayMass_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC("");
  for (size_t j(0);j<m_ids.size();++j) {
    Vec4D sum;
    for (size_t i(0);i<m_ids[j].size();++i) sum+=sl[m_ids[j][i]].Momentum();
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
    size_t id(ID(m_ids[j]));
    double scut(cuts->Getscut(id));
    cuts->Setscut(id,Max(scut,sqr(m_min)));
  }
}

DECLARE_GETTER(DecayMass_Selector,"DecayMass",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DecayMass_Selector>::
operator()(const Selector_Key &key) const
{
  DecayMass_Selector *msel(new DecayMass_Selector(key));
  return msel;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DecayMass_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[DecayMass, kf, min, max]";
}
