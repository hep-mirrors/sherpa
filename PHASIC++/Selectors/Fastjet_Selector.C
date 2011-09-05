#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "ATOOLS/Math/Algebra_Interpreter.H"

namespace PHASIC {
  class Fastjet_Selector: public Selector_Base, public ATOOLS::Tag_Replacer {
    double m_ptmin,m_etmin,m_delta_r,m_f;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;
    ATOOLS::Algebra_Interpreter m_calc;
    ATOOLS::Vec4D_Vector m_p;

    std::string ReplaceTags(std::string &expr) const;

    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  public:
    Fastjet_Selector(int nin, int nout,ATOOLS::Flavour * fl,std::string algo,
		     double ptmin, double etmin, double dr, double f,
		     int nn,std::string expression);

    ~Fastjet_Selector();


    bool   NoJetTrigger(const ATOOLS::Vec4D_Vector &);
    bool   Trigger(const ATOOLS::Vec4D_Vector &);
    bool   JetTrigger(const ATOOLS::Vec4D_Vector &,
		      ATOOLS::NLO_subevtlist *const subs);

    void   BuildCuts(Cut_Data *) {}
  };
}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"


using namespace PHASIC;
using namespace ATOOLS;

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Selector::Fastjet_Selector
(int nin, int nout,ATOOLS::Flavour * fl, std::string algo,
 double ptmin, double etmin, double dr, double f,
 int nn,std::string expression) : 
  Selector_Base("Fastjetfinder"), m_ptmin(ptmin), m_etmin(etmin), 
  m_delta_r(dr), m_f(f), p_jdef(0), p_siscplug(0)
{
  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);

  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(m_delta_r,m_f);

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r);

  m_fl         = fl;
  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));
  m_smax       = sqr(rpa->gen.Ecms());

  m_nin        = nin;
  m_nout       = nout;
  m_n          = nn;

  m_sel_log    = new Selector_Log(m_name);

  m_p.resize(m_n);

  m_calc.AddTag("H_T2","1.0");
  m_calc.AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  for (size_t i=0;i<m_p.size();++i) 
    m_calc.AddTag("p["+ToString(i)+"]",ToString(m_p[i]));

  m_calc.SetTagReplacer(this);
  m_calc.Interprete(expression);

  msg_Debugging()<<METHOD<<"(): '"<<expression<<"' {\n";
  msg_Indent();
  if (msg_LevelIsDebugging()) m_calc.PrintEquation();
  msg_Debugging()<<"}\n";
}


Fastjet_Selector::~Fastjet_Selector() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
}

std::string Fastjet_Selector::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Fastjet_Selector::ReplaceTags(Term *term) const
{
  if (term->Id()>=100) {
    term->Set(m_p[term->Id()-100]);
    return term;
  }
  else if (term->Id()==5) {
    double ht(0.0);
    for (size_t i(0);i<m_p.size();++i) ht+=m_p[i].PPerp();
    term->Set(sqr(ht));
    return term;
  }
  else if (term->Id()==6) {
    Vec4D sum(0.0,0.0,0.0,0.0);
    for (size_t i(0);i<m_p.size();++i) sum+=m_p[i];
    term->Set(sum);
    return term;
  }
  return term;
}

void Fastjet_Selector::AssignId(Term *term)
{
  if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag()=="P_SUM") term->SetId(6);
  else {
    int idx(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (idx>=m_n) THROW(fatal_error,"Invalid syntax");
    term->SetId(100+idx);
  }
}

bool Fastjet_Selector::NoJetTrigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  double s=(p[0]+p[1]).Abs2();
  return (s>m_smin*4.);
}

bool Fastjet_Selector::Trigger(const Vec4D_Vector &p)
{
  if (m_n<1) return true;

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<p.size();++i) {
    if (m_fl[i].Strong())
      input.push_back(fastjet::PseudoJet(p[i][1],p[i][2],p[i][3],p[i][0])); 
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  m_p.clear();
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin) m_p.push_back(pj);
  }

  bool trigger((int)m_p.size()>=m_n);
  if (trigger) trigger=(int)m_calc.Calculate()->Get<double>();

  return (1-m_sel_log->Hit(1-trigger));
}

bool Fastjet_Selector::JetTrigger(const Vec4D_Vector &p,
				ATOOLS::NLO_subevtlist *const subs)
{
  if (m_n<1) return true;

  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<subs->back()->m_n;++i) {
    if (subs->back()->p_fl[i].Strong())
      input.push_back(fastjet::PseudoJet(p[i][1],p[i][2],p[i][3],p[i][0]));      
  }
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=cs.inclusive_jets();

  m_p.clear();
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin) m_p.push_back(pj);
  }

  bool trigger((int)m_p.size()>=m_n);
  if (trigger) trigger=(int)m_calc.Calculate()->Get<double>();
  
  return (1-m_sel_log->Hit(1-trigger));
}


namespace PHASIC{

DECLARE_ND_GETTER(Fastjet_Selector_Getter,"FastjetSelector",Selector_Base,Selector_Key,true);

Selector_Base *Fastjet_Selector_Getter::operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<6) THROW(critical_error,"Invalid syntax");
 
  double f(.75);
  if (key.front().size()>6) f=ToType<double>(key[0][6]);

  Fastjet_Selector *jf(new Fastjet_Selector
		       (key.p_proc->NIn(),key.p_proc->NOut(),
			(Flavour*)&key.p_proc->Process()->Flavours().front(),key[0][1],
			ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3])),
			ToType<double>(key.p_read->Interpreter()->Interprete(key[0][4])),
			ToType<double>(key[0][5]),f,ToType<double>(key[0][2]),key[0][0]));
  jf->SetProcess(key.p_proc);
  return jf;
}

void Fastjet_Selector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Fastjet expression algorithm n ptmin etmin dr [f(siscone)=0.75]\n" 
     <<"        algorithm: kt,antikt,cambridge,siscone";
}

}

#endif
