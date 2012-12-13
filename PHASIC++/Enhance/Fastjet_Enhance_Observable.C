#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__FASTJET

#include "PHASIC++/Enhance/Enhance_Observable_Base.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/SISConePlugin.hh"

namespace PHASIC {

  class Fastjet_Enhance_Observable:
    public Enhance_Observable_Base,
    public ATOOLS::Tag_Replacer {
  private:

    ATOOLS::Algebra_Interpreter m_calc;

    ATOOLS::Vec4D_Vector m_p;
    std::vector<double> m_mu2;

    fastjet::JetDefinition *p_jdef;
    fastjet::SISConePlugin *p_siscplug;

    double m_ptmin, m_etmin, m_eta, m_y;

  public:

    Fastjet_Enhance_Observable(const Enhance_Arguments &args);

    double operator()(const ATOOLS::Vec4D *p,
		      const ATOOLS::Flavour *fl,
		      const size_t n);

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    void AssignId(ATOOLS::Term *term);

  };// end of class Fastjet_Enhance_Observable

}// end of namespace PHASIC

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Fastjet_Enhance_Observable_Getter,"FASTJET",
	       Enhance_Observable_Base,Enhance_Arguments);

Enhance_Observable_Base *Fastjet_Enhance_Observable_Getter::
operator()(const Enhance_Arguments &args) const
{
  return new Fastjet_Enhance_Observable(args);
}

void Fastjet_Enhance_Observable_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"fastjet enhance observable";
}

Fastjet_Enhance_Observable::Fastjet_Enhance_Observable
(const Enhance_Arguments &args): Enhance_Observable_Base(args),
  p_jdef(NULL), p_siscplug(NULL)
{
  std::string jtag(args.m_enhance);
  size_t pos(jtag.find("FASTJET["));
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid enhance function '"+args.m_enhance+"'");
  jtag=jtag.substr(pos+8);
  pos=jtag.find(']');
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid enhance function '"+args.m_enhance+"'");
  jtag=jtag.substr(0,pos);
  Data_Reader read(" ",",","#","=");
  read.AddIgnore(":");
  read.SetAddCommandLine(false);
  read.SetString(jtag);
  m_ptmin=read.StringValue<double>("PT",0.0);
  m_etmin=read.StringValue<double>("ET",0.0);
  m_eta=read.StringValue<double>("Eta",100.0);
  m_y=read.StringValue<double>("Y",100.0);
  double R(read.StringValue<double>("R",0.4));
  double f(read.StringValue<double>("f",0.75));
  std::string algo(read.StringValue<std::string>("A","antikt"));
  fastjet::JetAlgorithm ja(fastjet::kt_algorithm);
  if (algo=="cambridge") ja=fastjet::cambridge_algorithm;
  if (algo=="antikt") ja=fastjet::antikt_algorithm;
  if (algo=="siscone") p_siscplug=new fastjet::SISConePlugin(R,f);
  std::string reco(read.StringValue<std::string>("C","E"));
  fastjet::RecombinationScheme recom(fastjet::E_scheme);
  if (reco=="pt") recom=fastjet::pt_scheme;
  if (reco=="pt2") recom=fastjet::pt2_scheme;
  if (reco=="Et") recom=fastjet::Et_scheme;
  if (reco=="Et2") recom=fastjet::Et2_scheme;
  if (reco=="BIpt") recom=fastjet::BIpt_scheme;
  if (reco=="BIpt2") recom=fastjet::BIpt2_scheme;
  if (p_siscplug) p_jdef=new fastjet::JetDefinition(p_siscplug);
  else p_jdef=new fastjet::JetDefinition(ja,R,recom);
  m_mu2.resize(p_proc->NOut());
  m_p.resize(p_proc->NIn()+p_proc->NOut());
  std::string arg(args.m_enhance);
  size_t bpos(arg.find("{")), epos(arg.find("}",bpos));
  if (bpos==std::string::npos || epos==std::string::npos)
    THROW(fatal_error,"Invalid input");
  arg=arg.substr(bpos+1,arg.length()-bpos-2);
  m_calc.SetTagReplacer(this);
  for (int i(0);i<m_p.size();++i)
    m_calc.AddTag("p["+ToString(i)+"]",ToString(m_p[i]));
  for (size_t i=0;i<m_mu2.size();++i)
    m_calc.AddTag("MU_"+ToString(i)+"2",ToString(m_mu2[i]));
  m_calc.AddTag("H_T2","1.0");
  m_calc.Interprete(arg);
}

double Fastjet_Enhance_Observable::operator()
  (const ATOOLS::Vec4D *p,const ATOOLS::Flavour *fl,const size_t n)
{
  msg_Debugging()<<METHOD<<"("<<p_proc->Name()<<"): {\n";
  m_p.resize(0);
  for (size_t i(0);i<p_proc->NIn();++i) m_p.push_back(-p[i]);
  std::vector<fastjet::PseudoJet> input;
  for (size_t i(p_proc->NIn());i<n;++i) {
    msg_Debugging()<<"  p_"<<i<<" = "<<p[i]<<" ("<<fl[i]<<")\n";
    if (!fl[i].Strong()) m_p.push_back(p[i]);
    else input.push_back(fastjet::PseudoJet(p[i][1],p[i][2],p[i][3],p[i][0]));
  }
  fastjet::ClusterSequence cs(input,*p_jdef);
  std::vector<fastjet::PseudoJet>
    jets(fastjet::sorted_by_pt(cs.inclusive_jets()));
  for (size_t i(0);i<jets.size();++i) {
    Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
    if (pj.PPerp()>m_ptmin && pj.EPerp()>m_etmin &&
	(m_eta==100 || dabs(pj.Eta())<m_eta) &&
	(m_y==100 || dabs(pj.Y())<m_y)) {
      msg_Debugging()<<"  p_"<<m_p.size()<<" = "<<pj<<"\n";
      m_p.push_back(pj);
    }
  }
  m_mu2.resize(input.size());
  for (size_t i(0);i<input.size();++i) {
    m_mu2[i]=cs.exclusive_dmerge_max(i);
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_mu2[i])<<"\n";
  }
  double val(m_calc.Calculate()->Get<double>());
  msg_Debugging()<<"} -> "<<val<<"\n";
  return val;
}

std::string Fastjet_Enhance_Observable::ReplaceTags(std::string &expr) const
{
  return m_calc.ReplaceTags(expr);
}

Term *Fastjet_Enhance_Observable::ReplaceTags(Term *term) const
{
  if (term->Id()>=1000) {
    if (term->Id()-1000>=m_mu2.size())
      THROW(fatal_error,"\\mu index too large");
    term->Set(m_mu2[term->Id()-1000]);
    return term;
  }
  if (term->Id()>=100) {
    if (term->Id()-100>=m_p.size())
      THROW(fatal_error,"p index too large");
    term->Set(m_p[term->Id()-100]);
    return term;
  }
  if (term->Id()==5) {
    double ht(0.0);
    for (size_t i(0);i<m_p.size();++i) ht+=m_p[i].PPerp();
    term->Set(sqr(ht));
    return term;
  }
  return term;
}

void Fastjet_Enhance_Observable::AssignId(Term *term)
{
  if (term->Tag()=="H_T2") term->SetId(5);
  else if (term->Tag().find("MU_")==0) {
    int idx(ToType<int>(term->Tag().substr(3,term->Tag().length()-4)));
    if (idx>=m_mu2.size()) THROW(fatal_error,"Invalid syntax");
    term->SetId(1000+idx);
  }
  else {
    int idx(ToType<int>(term->Tag().substr(2,term->Tag().length()-3)));
    if (idx>=m_p.size()) THROW(fatal_error,"Invalid syntax");
    term->SetId(100+idx);
  }
}

#endif
