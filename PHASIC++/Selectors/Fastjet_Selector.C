#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "PHASIC++/Selectors/Selector.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/EECambridgePlugin.hh"
#include "fastjet/JadePlugin.hh"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"

namespace PHASIC {
  class Fastjet_Selector: public Selector_Base, public ATOOLS::Tag_Replacer {
    double m_ptmin,m_etmin,m_delta_r,m_f,m_eta,m_y;
    int m_nj, m_bmode, m_eekt;
    fastjet::JetDefinition * p_jdef;
    fastjet::SISConePlugin * p_siscplug;
    fastjet::EECambridgePlugin * p_eecamplug;
    fastjet::JadePlugin * p_jadeplug;
    ATOOLS::Algebra_Interpreter m_calc;
    ATOOLS::Vec4D_Vector m_p;
    std::vector<double> m_mu2;

    std::string ReplaceTags(std::string &expr) const;

    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  public:
    Fastjet_Selector(Process_Base *const proc, std::string algo, size_t nj,
                     double ptmin, double etmin, double dr, double f,
                     double eta, double y, int bmode,
                     std::string expression);

    ~Fastjet_Selector();


    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(Cut_Data *) {}
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

/*---------------------------------------------------------------------

  General form - flavours etc are unknown, will operate on a Particle_List.

  --------------------------------------------------------------------- */

Fastjet_Selector::Fastjet_Selector
(Process_Base *const proc, std::string algo, size_t nj,
 double ptmin, double etmin, double dr, double f, double eta, double y,
 int bmode,std::string expression) :
  Selector_Base("FastjetSelector",proc),
  m_nj(nj), m_ptmin(ptmin), m_etmin(etmin),
  m_delta_r(dr), m_f(f), m_eta(eta), m_y(y), m_bmode(bmode), m_eekt(0), p_jdef(0),
  p_siscplug(NULL), p_eecamplug(NULL), p_jadeplug(NULL)
{
  bool ee(rpa->gen.Beam1().IsLepton() && rpa->gen.Beam2().IsLepton());

  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);

  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt")    ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(m_delta_r,m_f);
  if (ee) {
    if (algo=="eecambridge") p_eecamplug=new fastjet::EECambridgePlugin(dr);
    if (algo=="jade") p_jadeplug=new fastjet::JadePlugin();
  }

  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else if (p_eecamplug) p_jdef=new fastjet::JetDefinition(p_eecamplug);
  else if (p_jadeplug) p_jdef=new fastjet::JetDefinition(p_jadeplug);
  else if (ee) {
    p_jdef=new fastjet::JetDefinition(fastjet::ee_kt_algorithm);
    m_eekt=1;
  }
  else p_jdef=new fastjet::JetDefinition(ja,m_delta_r);

  m_smin       = Max(sqr(m_ptmin),sqr(m_etmin));

  m_p.resize(m_nin+m_nout);
  m_mu2.resize(m_nout);

  m_calc.AddTag("H_T2","1.0");
  m_calc.AddTag("P_SUM","(1.0,0.0,0.0,0.0)");
  for (size_t i=0;i<m_p.size();++i)
    m_calc.AddTag("p["+ToString(i)+"]",ToString(m_p[i]));
  for (size_t i=0;i<m_mu2.size();++i)
    m_calc.AddTag("MU_"+ToString(i)+"2",ToString(m_mu2[i]));

  m_calc.SetTagReplacer(this);
  m_calc.Interprete(expression);

  msg_Debugging()<<METHOD<<"(): '"<<expression<<"' {\n";
  msg_Indent();
  if (msg_LevelIsDebugging()) m_calc.PrintEquation();
  msg_Debugging()<<"}\n";

#ifndef USING__FASTJET__3
  if (m_bmode>0) THROW(fatal_error, "b-tagging needs FastJet >= 3.0.");
#endif
}


Fastjet_Selector::~Fastjet_Selector() {
  delete p_jdef;
  if (p_siscplug) delete p_siscplug;
  if (p_eecamplug) delete p_eecamplug;
  if (p_jadeplug) delete p_jadeplug;
}

std::string Fastjet_Selector::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Fastjet_Selector::ReplaceTags(Term *term) const
{
  if (term->Id()>=1000) {
    term->Set(m_mu2[term->Id()-1000]);
    return term;
  }
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
  else if (term->Tag().find("MU_")==0) {
    int idx(ToType<int>(term->Tag().substr(3,term->Tag().length()-4)));
    if (idx>=m_mu2.size()) THROW(fatal_error,"Invalid syntax");
    term->SetId(1000+idx);
  }
  else {
    int idx(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (idx>=m_nin+m_nout) THROW(fatal_error,"Invalid syntax");
    term->SetId(100+idx);
  }
}

bool Fastjet_Selector::Trigger(Selector_List &sl)
{
  if (m_nj<0) return true;

  m_p.clear();
  for (size_t i(0);i<m_nin;++i) m_p.push_back(sl[i].Momentum());
  std::vector<fastjet::PseudoJet> input,jets;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (ToBeClustered(sl[i].Flavour(), m_bmode)) {
      input.push_back(MakePseudoJet(sl[i].Flavour(),sl[i].Momentum()));
    } else {
      m_p.push_back(sl[i].Momentum());
    }
  }
  int nj=m_p.size();
  
  fastjet::ClusterSequence cs(input,*p_jdef);
  jets=fastjet::sorted_by_pt(cs.inclusive_jets());

  if (m_eekt) {
    for (size_t i(0);i<input.size();++i) {
      if (cs.exclusive_dmerge_max(i)>sqr(m_ptmin)) {
        m_p.emplace_back(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
      }
    }
  } else {
    for (size_t i(0);i<jets.size();++i) {
      if (m_bmode==0 || BTag(jets[i], m_bmode)) {
        Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
        if (pj.PPerp()>m_ptmin&&pj.EPerp()>m_etmin &&
            (m_eta==100 || dabs(pj.Eta())<m_eta) &&
            (m_y==100 || dabs(pj.Y())<m_y)) m_p.push_back(pj);
      }
    }
  }
  for (size_t i(0);i<input.size();++i)
    m_mu2[i]=cs.exclusive_dmerge_max(i);

  bool trigger((int)(m_p.size()-nj)>=m_nj);
  if (trigger) trigger=(int)m_calc.Calculate()->Get<double>();

  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_ND_GETTER(Fastjet_Selector,"FastjetSelector",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Selector>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings;
  const auto expression = s["Expression"].SetDefault("").Get<std::string>();
  const auto algo = s["Algorithm"].SetDefault("").Get<std::string>();
  const auto n = s["N"].SetDefault(0).Get<size_t>();
  const auto ptmin = s["PTMin"].SetDefault(0.0).Get<double>();
  const auto etmin = s["ETMin"].SetDefault(0.0).Get<double>();
  const auto dr = s["DR"].SetDefault(0.4).Get<double>();
  const auto f = s["f"].SetDefault(0.75).Get<double>();
  const auto eta = s["EtaMax"].SetDefault(100.0).Get<double>();
  const auto y = s["YMax"].SetDefault(100.0).Get<double>();
  const auto bmode = s["BMode"].SetDefault(0).Get<int>();
  Fastjet_Selector *jf(
      new Fastjet_Selector(key.p_proc, algo, n, ptmin, etmin, dr,
                           f, eta, y, bmode, expression));
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Fastjet_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{\n"
     <<width<<"  Type: FastjetFinder,\n"
     <<width<<"  Expression: boolean expression,\n"
     <<width<<"  Algorithm: kt (default)|antikt|cambridge|siscone, (hadron colliders)\n"
     <<width<<"  Algorithm: eekt (default)|jade|eecambridge|siscone, (lepton colliders)\n"
     <<width<<"  N: number of jets,\n"
     <<width<<"  # optional settings:\n"
     <<width<<"  PTMin: minimum jet pT,\n"
     <<width<<"  ETMin: minimum jet eta,\n"
     <<width<<"  DR: jet distance parameter,\n"
     <<width<<"  f: Siscone f parameter,\n"
     <<width<<"  EtaMax: maximum jet eta,\n"
     <<width<<"  YMax: maximum jet rapidity,\n"
     <<width<<"  BMode: 0|1|2\n"
     <<width<<"  }";
}

#endif
