#include "PHASIC++/Selectors/Jet_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Jet_Finder::Jet_Finder
(Process_Integrator *const proc,
 const int nin,const int nout,Flavour *fl,
 const std::string &ycut):
  Selector_Base("Jetfinder"), m_cuttag(ycut),
  m_on(true), p_yccalc(NULL)
{
  p_proc=proc;
  m_ycut=2.0;
  m_fl=fl;
  m_nin=nin;
  m_nout=nout;
  m_n=m_nin+m_nout;
  m_smax=m_s=sqr(rpa->gen.Ecms());
  m_sel_log = new Selector_Log(m_name);
  static bool mets(false);
  if (!mets) {
    mets=true;
    rpa->gen.AddCitation(1,"LO/LO matrix element merging with truncated showers (MEPS/CKKW) is "+
			std::string("published under \\cite{Hoeche:2009rj}."));
  }
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetNIn(m_nin);
  for (int i(0);i<m_nin+m_nout;++i)
    p_ampl->CreateLeg(Vec4D(),i<m_nin?m_fl[i].Bar():m_fl[i],ColorID());
  p_ampl->SetJF(this);
  p_ampl->SetMS(proc->Process()->Generator());
  p_yccalc = new Algebra_Interpreter();
  p_yccalc->SetTagReplacer(this);
  for (int i=0;i<m_n;++i) p_yccalc->AddTag
    ("p["+ToString(i)+"]",ToString(Vec4D()));
  p_yccalc->Interprete(m_cuttag);
  p_jc = JetCriterion_Getter::GetObject
    (rpa->gen.Variable("JET_CRITERION"),
     JetCriterion_Key(rpa->gen.Variable("JET_CRITERION"),
		      p_proc->Process()->Shower()));
  if (p_jc==NULL) THROW(not_implemented,"Invalid jet criterion");
}

Jet_Finder::~Jet_Finder() 
{
  if (p_ampl) p_ampl->Delete();
  if (p_yccalc) delete p_yccalc;
  if (p_jc) delete p_jc;
}

bool Jet_Finder::Trigger(const Vec4D_Vector &p)
{
  p_ampl->SetProc(p_proc->Process());
  for (size_t i(0);i<p.size();++i)
    p_ampl->Leg(i)->SetMom((int)i<m_nin?-p[i]:p[i]);
  m_ycut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Process()->Name()
		 <<"' Q_cut = "<<sqrt(m_ycut*m_s)<<(m_on?" {":", off")<<"\n";
  p_ampl->Decays()=p_proc->Process()->Info().m_fi.GetDecayInfos();
  bool res=p_jc->Jets(p_ampl);
  msg_Debugging()<<"} -> "<<res<<"\n";
  return 1-m_sel_log->Hit(!res);
}

bool Jet_Finder::JetTrigger(const ATOOLS::Vec4D_Vector &p,
                            NLO_subevtlist *const subs)
{
  if (p.size()!=subs->back()->m_n) THROW(fatal_error,"Invalid call");
  p_ampl->SetProc(p_proc->Process());
  if (p_ampl->Legs().size()<subs->back()->m_n)
    p_ampl->CreateLeg(Vec4D(),Flavour(kf_jet),ColorID());
  else if (p_ampl->Legs().size()>subs->back()->m_n) {
    p_ampl->Legs().back()->Delete();
    p_ampl->Legs().pop_back();
  }
  size_t idij((1<<subs->back()->m_i)|(1<<subs->back()->m_j));
  if (subs->back()->m_i==subs->back()->m_j) idij=0;
  for (size_t i(0);i<p.size();++i) {
    p_ampl->Leg(i)->SetFlav
      ((int)i<m_nin?subs->back()->p_fl[i].Bar():subs->back()->p_fl[i]);
    p_ampl->Leg(i)->SetMom((int)i<m_nin?-p[i]:p[i]);
    p_ampl->Leg(i)->SetId(subs->back()->p_id[i]);
    p_ampl->Leg(i)->SetK(subs->back()->p_id[i]==idij?
			 (1<<subs->back()->m_k):0);
  }
  m_ycut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Process()->Name()
		 <<"' Q_cut = "<<sqrt(m_ycut*m_s)<<(m_on?" {":", off")<<"\n";
  p_ampl->Decays()=p_proc->Process()->Info().m_fi.GetDecayInfos();
  bool res=p_jc->Jets(p_ampl,idij?0:1);
  msg_Debugging()<<"} -> "<<res<<"\n";
  return 1-m_sel_log->Hit(!res);
}

bool Jet_Finder::NoJetTrigger(const ATOOLS::Vec4D_Vector &p)
{
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

DECLARE_ND_GETTER(Jet_Finder,"METS",Selector_Base,Selector_Key,false);
  
Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Jet_Finder>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<1) THROW(critical_error,"Invalid syntax");
  Jet_Finder *jf(new Jet_Finder(key.p_proc,key.p_proc->NIn(),key.p_proc->NOut(),
				(Flavour*)&key.p_proc->Process()->
				Flavours().front(),key[0][0]));
  static bool menlots(false);
  if (!menlots && key.p_proc->Process()->Info().Has(nlo_type::real)) {
    menlots=true;
    rpa->gen.AddCitation(1,"NLO/LO matrix element merging with truncated showers (MENLOPS) is "+
			 std::string("published under \\cite{Hoeche:2010kg}."));
    rpa->gen.AddCitation(1,"NLO/NLO matrix element merging with truncated showers (MEPS@NLO) is "+
                         std::string("published under \\cite{Hoeche:2012yf} and \\cite{Gehrmann:2012yg}."));
  }
  if (key.front().size()>1 && key[0][1]=="LO" && 
      !(key.front().size()>2 && key[0][2]=="CUT")) 
    jf->SetOn(false);
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Jet_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"METS jet finder"; 
}

