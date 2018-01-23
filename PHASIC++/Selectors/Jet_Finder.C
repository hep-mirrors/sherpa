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
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"

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
  p_jc = JetCriterion_Getter::GetObject
    (rpa->gen.Variable("JET_CRITERION"),
     JetCriterion_Key(rpa->gen.Variable("JET_CRITERION"),p_proc->Shower()));
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
  p_ampl->SetProc(p_proc);
  for (size_t i(0);i<sl.size();++i)
    p_ampl->Leg(i)->SetMom((int)i<m_nin?-sl[i].Momentum():sl[i].Momentum());
  m_qcut=p_yccalc->Calculate()->Get<double>();
  if (!m_on) return true;
  msg_Debugging()<<METHOD<<"("<<this<<"): '"<<p_proc->Name()
		 <<"' Q_cut = "<<m_qcut<<(m_on?" {":", off")<<"\n";
  p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
  double jcv=p_jc->Value(p_ampl);
  bool res=m_pass=jcv>sqr(m_qcut);
  msg_Debugging()<<"} -> "<<res<<"\n";
  if (p_proc->VariationWeights()) {
    Reweight_Args args(jcv,0);
    p_proc->VariationWeights()->UpdateOrInitialiseWeights
      (&Jet_Finder::Reweight,*this,args);
    if (args.m_acc) res=1;
  }
  return 1-m_sel_log->Hit(!res);
}

double Jet_Finder::Reweight(Variation_Parameters *params,
			    Variation_Weights *weights,
			    Reweight_Args &args)
{
  msg_Debugging()<<METHOD<<"(): '"<<p_proc->Name()
		 <<"' Q_cut = "<<m_qcut*params->m_Qcutfac<<"\n";
  bool res=args.m_jcv>sqr(m_qcut*params->m_Qcutfac);
  msg_Debugging()<<"  jcv = "<<sqrt(args.m_jcv)<<"\n";
  msg_Debugging()<<"} -> "<<res<<"\n";
  if (res) args.m_acc=1;
  return res?1.0:0.0;
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
  Jet_Finder *jf(new Jet_Finder(key.p_proc,key[0][0]));
  static bool menlots(false);
  if (!menlots && key.p_proc->Info().Has(nlo_type::real)) {
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

