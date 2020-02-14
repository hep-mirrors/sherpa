#include "PHASIC++/Selectors/Jet_Finder.H"

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Jet_Finder::Jet_Finder(Process_Base *const proc,const std::string &ycut):
  Selector_Base("Jetfinder",proc), m_cuttag(ycut),
  m_on(true), p_yccalc(NULL)
{
  static bool mets(false);
  if (!mets) {
    mets=true;
    rpa->gen.AddCitation(1,"LO/LO matrix element merging with truncated showers (MEPS/CKKW) is "+
			std::string("published under \\cite{Hoeche:2009rj}."));
  }
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetNIn(m_nin);
  for (int i(0);i<m_nin+m_nout;++i)
    p_ampl->CreateLeg(Vec4D(),i<m_nin?p_fl[i].Bar():p_fl[i],ColorID());
  p_ampl->SetJF(this);
  p_ampl->SetMS(proc->Generator());
  p_yccalc = new Algebra_Interpreter();
  p_yccalc->SetTagReplacer(this);
  for (int i=0;i<m_n;++i) p_yccalc->AddTag
    ("p["+ToString(i)+"]",ToString(Vec4D()));
  p_yccalc->Interprete(m_cuttag);
  Settings& s = Settings::GetMainSettings();
  p_jc = JetCriterion_Getter::GetObject(
      s["JET_CRITERION"].Get<std::string>(),
      JetCriterion_Key(s["JET_CRITERION"].Get<std::string>(),
                       p_proc->Shower()));
  if (p_jc==NULL) THROW(not_implemented,"Invalid jet criterion");
}

Jet_Finder::~Jet_Finder() 
{
  if (p_ampl) p_ampl->Delete();
  if (p_yccalc) delete p_yccalc;
  if (p_jc) delete p_jc;
}

bool Jet_Finder::Trigger(Selector_List &sl)
{
  m_pass=true;
  m_results = {{1.0}};
  p_ampl->SetProc(p_proc);
  for (size_t i(0);i<sl.size();++i)
    p_ampl->Leg(i)->SetMom((int)i<m_nin?-sl[i].Momentum():sl[i].Momentum());
  m_qcut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"("<<this<<"): '"<<p_proc->Name()
		 <<"' Q_cut = "<<m_qcut<<(m_on?" {":", off")<<"\n";
  p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
  const double jcv=p_jc->Value(p_ampl);
  bool triggered {false};
  m_results[0].ApplyAll(
      [this, jcv, &triggered](double varweight,
                              size_t varindex,
                              Variation_Parameters* varparams) -> double {
        const bool pass =
            jcv > sqr(m_qcut * (varparams ? varparams->m_Qcutfac : 1.0));
        if (!varparams)
          m_pass = pass;
        if (pass)
          triggered = true;
        return pass;
      });
  msg_Debugging()<<"} -> "<<m_pass<<"\n";
  return 1-m_sel_log->Hit(!triggered);
}

bool Jet_Finder::RSTrigger(NLO_subevtlist *const subs)
{
  for (size_t i(0);i<m_nin+m_nout;++i)
    p_ampl->Leg(i)->SetMom(i<m_nin && subs->back()->p_mom[i][0]>0.0?
			   -subs->back()->p_mom[i]:subs->back()->p_mom[i]);
  m_qcut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  m_pass=0;
  std::vector<int> any_variation_passes(subs->size(), 0);
  m_results = std::vector<Event_Weights>(subs->size(), {0.0});
  for (size_t n(0);n<subs->size();++n) {
    int nominal_passes {0};
    msg_Debugging()<<METHOD<<"("<<n<<"): '"<<p_proc->Name()
		   <<"' Q_cut = "<<m_qcut<<(m_on?" {":", off")<<"\n";
    {
      msg_Indent();
      p_ampl->SetProc(p_proc);
      if (p_ampl->Legs().size()<(*subs)[n]->m_n)
	p_ampl->CreateLeg(Vec4D(),Flavour(kf_jet),ColorID());
      else if (p_ampl->Legs().size()>(*subs)[n]->m_n) {
	p_ampl->Legs().back()->Delete();
	p_ampl->Legs().pop_back();
      }
      size_t idij((1<<(*subs)[n]->m_i)|(1<<(*subs)[n]->m_j));
      if ((*subs)[n]->m_i==(*subs)[n]->m_j) idij=0;
      for (size_t i(0);i<(*subs)[n]->m_n;++i) {
	p_ampl->Leg(i)->SetFlav
	  ((int)i<m_nin?(*subs)[n]->p_fl[i].Bar():(*subs)[n]->p_fl[i]);
	p_ampl->Leg(i)->SetMom(i<m_nin && subs->back()->p_mom[i][0]>0.0?
			       -(*subs)[n]->p_mom[i]:(*subs)[n]->p_mom[i]);
	p_ampl->Leg(i)->SetId((*subs)[n]->p_id[i]);
	p_ampl->Leg(i)->SetK((*subs)[n]->p_id[i]==idij?
			     (1<<(*subs)[n]->m_k):0);
      }
      p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
      const double jcv = p_jc->Value(p_ampl, idij ? 0 : 1);
      m_results[n].ApplyAll(
          [this, &any_variation_passes, &nominal_passes, &n, subs, jcv](
              double varweight,
              size_t varindex,
              Variation_Parameters* varparams) -> double {
            const auto pass =
                jcv > sqr(m_qcut * (varparams ? varparams->m_Qcutfac : 1.0));
            if (!varparams) {
              (*subs)[n]->m_trig = pass;
              if (pass)
                nominal_passes = m_pass = pass;
            }
            if (pass) {
              any_variation_passes[n] = 1;
            }
            return pass;
          });
    }
    msg_Debugging()<<"} -> "<<nominal_passes<<"\n";
  }

  int result = m_pass;
  for (size_t n(0); n < subs->size(); ++n) {
    result |= (*subs)[n]->m_trig |= (any_variation_passes[n] ? 2 : 0);
  }

  return 1-m_sel_log->Hit(!result);
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
  auto s = key.m_settings["METS"];
  const auto ycut = s["YCUT"].SetDefault("").Get<std::string>();
  if (ycut == "")
    THROW(critical_error,"Invalid syntax");
  auto* jf = new Jet_Finder(key.p_proc, ycut);
  static bool menlots(false);
  if (!menlots && key.p_proc->Info().Has(nlo_type::real)) {
    menlots=true;
    rpa->gen.AddCitation(1,"NLO/LO matrix element merging with truncated showers (MENLOPS) is "+
			 std::string("published under \\cite{Hoeche:2010kg}."));
    rpa->gen.AddCitation(1,"NLO/NLO matrix element merging with truncated showers (MEPS@NLO) is "+
                         std::string("published under \\cite{Hoeche:2012yf} and \\cite{Gehrmann:2012yg}."));
  }
  if (s["LO"].SetDefault(false).Get<bool>()
      && !s["CUT"].SetDefault(false).Get<bool>()) {
    jf->SetOn(false);
  }
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Jet_Finder>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"METS jet finder"; 
}
