#include "PHASIC++/Selectors/Jet_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace ATOOLS;

Jet_Finder::Jet_Finder
(const int nin,const int nout,Flavour *fl,
 const std::string &ycut):
  Selector_Base("Jetfinder"), m_dparam(0.3), m_cuttag(ycut),
  m_on(true), p_yccalc(NULL)
{
  m_ycut=2.0;
  m_fl=fl;
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_smax=m_s=sqr(rpa.gen.Ecms());
  if (ycut.find("|")!=std::string::npos) {
    m_dparam=ToType<double>(ycut.substr(ycut.find("|")+1));
    m_cuttag=ycut.substr(0, ycut.find("|"));
  }
  m_sel_log = new Selector_Log(m_name);
  static bool mets(false);
  if (!mets) {
    mets=true;
    rpa.gen.AddCitation(1,"Matrix element merging with truncated showers is "+
			std::string("published under \\cite{Hoeche:2009rj}."));
  }
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetNIn(m_nin);
  for (int i(0);i<m_nin+m_nout;++i)
    p_ampl->CreateLeg(Vec4D(),i<m_nin?m_fl[i].Bar():m_fl[i],ColorID());
  p_ampl->SetJF(this);
  p_yccalc = new Algebra_Interpreter();
  p_yccalc->SetTagReplacer(this);
  for (int i=0;i<m_n;++i) p_yccalc->AddTag
    ("p["+ToString(i)+"]",ToString(Vec4D()));
  p_yccalc->Interprete(m_cuttag);
}

Jet_Finder::~Jet_Finder() 
{
  p_ampl->Delete();
  delete p_yccalc;
}

bool Jet_Finder::Trigger(const Vec4D_Vector &p)
{
  for (size_t i(0);i<p.size();++i)
    p_ampl->Leg(i)->SetMom((int)i<m_nin?-p[i]:p[i]);
  m_ycut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Process()->Name()
		 <<"' Q_cut = "<<sqrt(m_ycut*m_s)<<(m_on?" {":", off")<<"\n";
  bool res=p_proc->Process()->Shower()->JetVeto(p_ampl);
  msg_Debugging()<<"} -> "<<res<<"\n";
  return 1-m_sel_log->Hit(!res);
}

bool Jet_Finder::JetTrigger(const ATOOLS::Vec4D_Vector &p,
                            NLO_subevtlist *const subs)
{
  if (!m_on) return true;
  THROW(not_implemented,"Don't even try!");
  return true;
}

bool Jet_Finder::NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
{
  if (!m_on) return true;
  THROW(not_implemented,"Don't even try!");
  return true;
}

void Jet_Finder::BuildCuts(Cut_Data *cuts) 
{
}

std::string Jet_Finder::ReplaceTags(std::string &expr) const
{
  return p_yccalc->ReplaceTags(expr);
}

Term *Jet_Finder::ReplaceTags(Term *term) const
{
  term->Set(p_ampl->Leg(term->Id())->Mom());
  return term;
}

void Jet_Finder::AssignId(Term *term)
{
  term->SetId(ToType<int>
	      (term->Tag().substr
	       (2,term->Tag().length()-3)));
}

namespace PHASIC{
  
  DECLARE_ND_GETTER(Jet_Finder_Getter,"METS",Selector_Base,Selector_Key,false);
  
  Selector_Base *Jet_Finder_Getter::operator()(const Selector_Key &key) const
  {
    if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
    Jet_Finder *jf(new Jet_Finder(key.p_proc->NIn(),key.p_proc->NOut(),
				  (Flavour*)&key.p_proc->Process()->
				  Flavours().front(),key[0][0]));
    jf->SetProcess(key.p_proc);
    static bool menlots(false);
    if (!menlots && key.p_proc->Process()->Info().Has(nlo_type::real)) {
      menlots=true;
      rpa.gen.AddCitation(1,"NLO matrix element merging with truncated showers is "+
			  std::string("published under \\cite{Hoeche:2010kg}."));
    }
    if (key.front().size()>1 && key[0][1]=="LO" && 
	!(key.front().size()>2 && key[0][2]=="CUT")) 
      jf->SetOn(false);
    return jf;
  }
  
  void Jet_Finder_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"METS jet finder"; 
  }
  
}
