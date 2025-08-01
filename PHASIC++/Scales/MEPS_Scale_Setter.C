#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PHASIC++/Scales/Color_Setter.H"
#include "PHASIC++/Scales/Cluster_Definitions.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

namespace PHASIC {

  class MEPS_Scale_Setter: public Scale_Setter_Base {
  private:

    PDF::Cluster_Definitions_Base *p_clu, *p_qdc;
    Core_Scale_Setter *p_core, *p_uoscale;
    Color_Setter *p_cs;
    ATOOLS::Mass_Selector *p_ms;
    Single_Process *p_sproc;

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

    std::shared_ptr<Color_Integrator> p_ci;

    double m_rsf, m_fsf;
    int    m_cmode, m_nmin;
    int    m_rproc, m_sproc, m_rsproc, m_vproc, m_nproc;

    static int s_nfgsplit, s_allowuo, s_nlocpl;

    int Select(const PDF::ClusterInfo_Vector &ccs,
	       const Int_Vector &on,const int mode=0) const;

    bool CoreCandidate(ATOOLS::Cluster_Amplitude *const ampl) const;
    bool CheckOrdering(ATOOLS::Cluster_Amplitude *const ampl,
		       const int ord=1) const;
    bool CheckSplitting(const PDF::Cluster_Info &cp,
			const int ord) const;
    bool CheckSubEvents(const PDF::Cluster_Config &cc) const;
    void Combine(ATOOLS::Cluster_Amplitude &ampl,
		 const PDF::Cluster_Info &ci) const;

    double SetScales(ATOOLS::Cluster_Amplitude *ampl);
    double Differential(ATOOLS::Cluster_Amplitude *const ampl,
			const int mode=0) const;

    bool ClusterStep(ATOOLS::Cluster_Amplitude *ampl,
		     ATOOLS::ClusterAmplitude_Vector &ampls,
		     const PDF::Cluster_Info &ci,const int ord) const;
    void Cluster(ATOOLS::Cluster_Amplitude *ampl,
		 ATOOLS::ClusterAmplitude_Vector &ampls,
		 const int ord) const;
    void SetCoreScale(ATOOLS::Cluster_Amplitude *const ampl) const;

    double UnorderedScale(ATOOLS::Cluster_Amplitude *const ampl) const;

  public:

    MEPS_Scale_Setter(const Scale_Setter_Arguments &args,
		      const int mode=1);

    ~MEPS_Scale_Setter();

    double Calculate(const ATOOLS::Vec4D_Vector &p,const size_t &mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

  };// end of class MEPS_Scale_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

DECLARE_GETTER(MEPS_Scale_Setter,"METS",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,MEPS_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new MEPS_Scale_Setter(args,1);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    MEPS_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"meps scale scheme";
}

int MEPS_Scale_Setter::s_nfgsplit(-1);
int MEPS_Scale_Setter::s_allowuo(0);
int MEPS_Scale_Setter::s_nlocpl(-1);

MEPS_Scale_Setter::MEPS_Scale_Setter
(const Scale_Setter_Arguments &args,const int mode):
  Scale_Setter_Base(args), m_tagset(this)
{
  static std::string s_core, s_uoscale;
  static int s_cmode(-1), s_csmode, s_nmaxall, s_nmaxnloall, s_kfac;
  if (s_cmode<0) {
    Scoped_Settings s(Settings::GetMainSettings()["MEPS"]);
    s_nmaxall=s["NMAX_ALLCONFIGS"].GetScalarWithOtherDefault<int>(-1);
    s_nmaxnloall=s["NLO_NMAX_ALLCONFIGS"].GetScalarWithOtherDefault<int>(-1);
    s_cmode=s["CLUSTER_MODE"].GetScalarWithOtherDefault<int>(8|32|64);
    s_nlocpl=s["NLO_COUPLING_MODE"].GetScalarWithOtherDefault<int>(2);
    s_csmode=s["MEPS_COLORSET_MODE"].GetScalarWithOtherDefault<int>(0);
    s_core=s["CORE_SCALE"].GetScalarWithOtherDefault<std::string>("Default");
    s_allowuo=s["ALLOW_SCALE_UNORDERING"].GetScalarWithOtherDefault<int>(0);
    s_uoscale=s["UNORDERED_SCALE"].GetScalarWithOtherDefault<std::string>("None");
    s_nfgsplit=Settings::GetMainSettings()["DIPOLES"]["NF_GSPLIT"].Get<int>();
    s_kfac = Settings::GetMainSettings()["SHOWER"]["KFACTOR_SCHEME"].Get<int>();
  }
  m_scale.resize(2*stp::size);
  std::string tag(args.m_scale), core(s_core), uoscale(s_uoscale);
  m_nproc=!(p_proc->Info().m_fi.NLOType()==nlo_type::lo);
  m_nmin=p_proc->Info().m_fi.NMinExternal();
  size_t pos(tag.find('['));
  if (pos==4) {
    tag=tag.substr(pos+1);
    pos=tag.find(']');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    Data_Reader read(" ",",","#","=");
    read.AddIgnore(":");
    read.SetString(tag.substr(0,pos));
    core=read.StringValue<std::string>("C",core);
    uoscale=read.StringValue<std::string>("U",uoscale);
    m_nmin=read.StringValue<int>("M",m_nmin);
    tag=tag.substr(pos+1);
  }
  if (tag.find('{')==std::string::npos &&
      tag.find('}')==std::string::npos) tag+="{MU_F2}{MU_R2}{MU_Q2}";
  while (true) {
    size_t pos(tag.find('{'));
    if (pos==std::string::npos) {
      if (!m_calcs.empty()) break;
      else { THROW(fatal_error,"Invalid scale '"+args.m_scale+"'"); }
    }
    tag=tag.substr(pos+1);
    pos=tag.find('}');
    if (pos==std::string::npos)
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    std::string ctag(tag.substr(0,pos));
    tag=tag.substr(pos+1);
    pos=ctag.find('|');
    if (pos!=std::string::npos)
      ctag=m_nproc?ctag.substr(0,pos):ctag.substr(pos+1);
    m_calcs.push_back(new Algebra_Interpreter());
    m_calcs.back()->AddFunction(MODEL::as->GetAIGMeanFunction());
    m_calcs.back()->SetTagReplacer(&m_tagset);
    if (m_calcs.size()==1) m_tagset.SetCalculator(m_calcs.back());
    SetScale(ctag,*m_calcs.back());
  }
  if (p_proc->Shower()==NULL) {
    msg_Error()<<"\nPlease specify the renormalisation and factorisation scale e.g SCALES: VAR{H_T2}"<<endl
               <<"For detailed instructions on setting scales, please refer to the Sherpa manual."<<endl;
    THROW(missing_input,"No shower generator found. Unable to use the MEPS scale scheme.");
  }
  p_clu=p_proc->Shower()->GetClusterDefinitions();
  p_ms=p_proc->Generator();
  m_scale.resize(Max(m_scale.size(),m_calcs.size()));
  SetCouplings();
  int ne(p_proc->NOut()-p_proc->Info().m_fi.NMinExternal());
  int ac(ne<=s_nmaxall?2:0);
  if (m_nproc && ne<=s_nmaxnloall) ac=2;
  m_cmode=s_cmode|ac;
  /*
    1 - Approximate probabilities
    2 - Construct all configurations
    4 - Winner takes it all
    8 - Ignore color
    16 - Do not include incomplete paths
    32 - Winner takes it all at NLO
    64 - Winner takes it all in R first step
    128 - Use R configuration in all RS
    256 - No ordering check if last qcd split
    512 - No ordering check if first RS split
    1024 - No differential for core
    2048 - No tracking of scale sums
  */
  p_core=Core_Scale_Getter::GetObject(core,Core_Scale_Arguments(p_proc,core));
  if (p_core==NULL) THROW(fatal_error,"Invalid core scale '"+core+"'");
  if (s_uoscale=="None") p_uoscale=NULL;
  else {
    p_uoscale=Core_Scale_Getter::GetObject(uoscale,Core_Scale_Arguments(p_proc,uoscale));
    if (p_uoscale==NULL) THROW(fatal_error,"Invalid unordered scale '"+uoscale+"'");
    msg_Debugging()<<METHOD<<"(): Custom scale for unordered configurations '"+uoscale+"'\n";
  }
  m_rsf=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  if (m_rsf!=1.0)
    msg_Debugging()<<METHOD<<"(): Renormalization scale factor "<<sqrt(m_rsf)<<"\n";
  m_fsf=ToType<double>(rpa->gen.Variable("FACTORIZATION_SCALE_FACTOR"));
  if (m_fsf!=1.0)
    msg_Debugging()<<METHOD<<"(): Factorization scale factor "<<sqrt(m_fsf)<<"\n";
  p_cs = new Color_Setter(s_csmode);
  p_qdc = new Cluster_Definitions(s_kfac,m_nproc,p_proc->Shower()->KTType());
}

MEPS_Scale_Setter::~MEPS_Scale_Setter()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
  for (size_t i(0);i<m_ampls.size();++i) m_ampls[i]->Delete();
  if (p_uoscale) delete p_uoscale;
  delete p_core;
  delete p_cs;
  delete p_qdc;
}

int MEPS_Scale_Setter::Select
(const ClusterInfo_Vector &ccs,const Int_Vector &on,const int mode) const
{
  if (mode==1 || (m_cmode&4) || ((m_cmode&32) && m_nproc)) {
    int imax(-1);
    double max(0.0);
    for (size_t i(0);i<ccs.size();++i)
      if (on[i] && dabs(ccs[i].second.m_op)>max) {
	max=dabs(ccs[i].second.m_op);
	imax=i;
      }
    return imax;
  }
  double sum(0.0);
  for (size_t i(0);i<ccs.size();++i)
    if (on[i]) sum+=dabs(ccs[i].second.m_op);
  if (sum==0.0) return -1;
  double disc(sum*ran->Get());
  sum=0.0;
  for (size_t i(0);i<ccs.size();++i) if (on[i])
      if ((sum+=dabs(ccs[i].second.m_op))>=disc) return i;
  return -1;
}

bool MEPS_Scale_Setter::CheckOrdering
(Cluster_Amplitude *const ampl,const int ord) const
{
  if (ampl->Prev()==NULL) return true;
  if ((m_cmode&512) && m_rproc &&
      ampl->Prev()->Prev()==NULL) return true;
  if (ampl->KT2()<ampl->Prev()->KT2()) {
    if ((m_cmode&256) &&
	(ampl->OrderQCD()==(m_vproc?1:0) ||
	 (ampl->OrderQCD()>(m_vproc?2:1) && ampl->Legs().size()==3))) {
      msg_Debugging()<<"No ordering veto: "<<sqrt(ampl->KT2())
		     <<" < "<<sqrt(ampl->Prev()->KT2())<<"\n";
      return true;
    }
    msg_Debugging()<<"Veto ordering: "<<sqrt(ampl->KT2())
		   <<" < "<<sqrt(ampl->Prev()->KT2())<<"\n";
    return false;
  }
  return true;
}

bool MEPS_Scale_Setter::CheckSplitting
(const Cluster_Info &ci,const int ord) const
{
  if (!CheckOrdering(ci.first.p_ampl,ord)) return false;
  Cluster_Leg *li(ci.first.p_ampl->Leg(ci.first.m_i));
  Cluster_Leg *lj(ci.first.p_ampl->Leg(ci.first.m_j));
  if (ci.first.m_mo.IsGluon() &&
      !li->Flav().IsGluon() && li->Flav().Kfcode()>s_nfgsplit &&
      !lj->Flav().IsGluon() && lj->Flav().Kfcode()>s_nfgsplit) {
    msg_Debugging()<<"Veto flavour\n";
    return false;
  }
  if (ci.first.p_ampl->OrderQCD()<
      (ci.second.m_cpl?(ci.second.m_cpl&2):1) ||
      ci.first.p_ampl->OrderEW()<(ci.second.m_cpl?1:0)) {
    msg_Debugging()<<"Veto order\n";
    return false;
  }
  return true;
}

bool MEPS_Scale_Setter::CheckSubEvents(const Cluster_Config &cc) const
{
  NLO_subevtlist *subs(p_proc->Caller()->GetRSSubevtList());
  for (size_t i(0);i<subs->size()-1;++i) {
    NLO_subevt *sub((*subs)[i]);
    Flavour mofl(sub->p_fl[sub->m_ijt]);
    if (sub->m_ijt<p_proc->NIn()) mofl=mofl.Bar();
    if (cc.m_k==sub->m_k && cc.m_mo==mofl &&
	((cc.m_i==sub->m_i && cc.m_j==sub->m_j) ||
	 (cc.m_i==sub->m_j && cc.m_j==sub->m_i))) return true;
  }
  return false;
}

void MEPS_Scale_Setter::Combine
(Cluster_Amplitude &ampl,const Cluster_Info &ci) const
{
  int i(ci.first.m_i), j(ci.first.m_j);
  if (i>j) std::swap<int>(i,j);
  Cluster_Leg *li(ampl.Leg(i)), *lj(ampl.Leg(j));
  Cluster_Leg *lk(ampl.Leg(ci.first.m_k));
  li->SetCol(ampl.CombineColors(li,lj,lk,ci.first.m_mo));
  li->SetFlav(ci.first.m_mo);
  li->SetMom(ci.second.m_pijt);
  li->SetStat(ci.second.m_stat);
  if (!(m_cmode&2048)) {
    li->SetKT2(0,sqr(sqrt(ci.second.m_kt2)+
		     sqrt(li->KT2(0))+sqrt(lj->KT2(0))));
    li->SetKT2(1,ci.second.m_kt2+li->KT2(1)+lj->KT2(1));
  }
  lk->SetMom(ci.second.m_pkt);
  ampl.Prev()->SetIdNew(ampl.Leg(ci.first.m_j)->Id());
  for (size_t m(0);m<ampl.Legs().size();++m) {
    ampl.Leg(m)->SetK(0);
    ampl.Leg(m)->SetStat(ampl.Leg(m)->Stat()|1);
  }
  if (i<2) {
    for (size_t m(0);m<ampl.Legs().size();++m) {
      if (ampl.Prev()) ampl.Leg(m)->SetNMax(ampl.Prev()->Leg(m)->NMax());
      if ((int)m==i || (int)m==j || (int)m==ci.first.m_k) continue;
      ampl.Leg(m)->SetMom(ci.second.m_lam*ampl.Leg(m)->Mom());
    }
  }
  li->SetId(li->Id()+lj->Id());
  li->SetK(lk->Id());
  std::vector<Cluster_Leg*>::iterator lit(ampl.Legs().begin());
  for (int l(0);l<j;++l) ++lit;
  (*lit)->Delete();
  ampl.Legs().erase(lit);
  ampl.SetOrderQCD(ampl.OrderQCD()-(ci.second.m_cpl?0:1));
  if (ci.second.m_cpl==2) ampl.SetOrderQCD(ampl.OrderQCD()-2);
  ampl.SetOrderEW(ampl.OrderEW()-(ci.second.m_cpl?1:0));
  ampl.SetKin(ci.second.m_kin);
}

double MEPS_Scale_Setter::Calculate
(const Vec4D_Vector &momenta,const size_t &mode)
{
  m_p=momenta;
  if (m_nproc || (m_cmode&8)) p_ci=NULL;
  else p_ci=p_proc->Caller()->Integrator()->ColorIntegrator();
  for (size_t i(0);i<p_proc->Caller()->NIn();++i) m_p[i]=-m_p[i];
  while (m_ampls.size()) {
    m_ampls.back()->Delete();
    m_ampls.pop_back();
  }
  DEBUG_FUNC(p_proc->Name()<<" from "<<p_proc->Caller()->Name());
  m_rproc=p_proc->Caller()->Info().Has(nlo_type::real);
  m_sproc=p_proc->Caller()->Info().Has(nlo_type::rsub);
  m_vproc=p_proc->Info().Has(nlo_type::vsub);
  m_rsproc=m_rproc||(m_sproc&&p_proc->Caller()->GetRSSubevtList());
  msg_Debugging()<<"vproc = "<<m_vproc<<", sproc = "<<m_sproc
		 <<", rproc = "<<m_rproc<<", rsproc = "<<m_rsproc<<"\n";
  const Flavour_Vector &fl=p_proc->Caller()->Flavours();
  if ((m_cmode&128) && m_sproc && !m_rproc &&
      p_proc->Caller()->GetRSSubevtList()) {
    msg_Debugging()<<"S event, use R configuration\n";
    Cluster_Amplitude *ampl(p_proc->Caller()->GetRSSubevtList()->back()->p_ampl);
    if (ampl && ampl->Next()) {
      m_ampls.push_back(ampl->CopyAll());
      m_ampls.back()=m_ampls.back()->Next();
      m_ampls.back()->DeletePrev();
    }
    else {
      Cluster_Amplitude *ampl(Cluster_Amplitude::New());
      ampl->SetNIn(p_proc->Caller()->NIn());
      ampl->SetOrderQCD(p_proc->Caller()->MaxOrder(0));
      for (size_t i(1);i<p_proc->Caller()->MaxOrders().size();++i)
	ampl->SetOrderEW(ampl->OrderEW()+p_proc->Caller()->MaxOrder(i));
      for (size_t i(0);i<m_p.size();++i)
	ampl->CreateLeg(m_p[i],i<p_proc->NIn()?fl[i].Bar():fl[i]);
      ampl->SetProc(p_proc->Caller()->Get<Single_Process>());
      SetCoreScale(ampl);
      m_ampls.push_back(ampl);
    }
    return SetScales(m_ampls.back());
  }
  size_t nmax(p_proc->Info().m_fi.NMaxExternal());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  ampl->SetNIn(p_proc->Caller()->NIn());
  ampl->SetOrderQCD(p_proc->Caller()->MaxOrder(0));
  for (size_t i(1);i<p_proc->Caller()->MaxOrders().size();++i)
    ampl->SetOrderEW(ampl->OrderEW()+p_proc->Caller()->MaxOrder(i));
  ampl->SetJF(p_proc->Selector()->GetSelector("Jetfinder"));
  if (p_ci != nullptr) {
    Int_Vector ci(p_ci->I()), cj(p_ci->J());
    for (size_t i(0);i<m_p.size();++i) {
      ampl->CreateLeg(m_p[i],i<p_proc->NIn()?fl[i].Bar():fl[i],
		      ColorID(ci[i],cj[i]));
      ampl->Leg(i)->SetNMax(nmax);
      ampl->Leg(i)->SetKT2(0,0.);
      ampl->Leg(i)->SetKT2(1,0.);
    }
  }
  else {
    for (size_t i(0);i<m_p.size();++i) {
      ampl->CreateLeg(m_p[i],i<p_proc->NIn()?fl[i].Bar():fl[i]);
      int cr(ampl->Leg(i)->Flav().StrongCharge());
      ampl->Leg(i)->SetCol(ColorID((cr==3||cr==8)?1:0,(cr==-3||cr==8)?1:0));
      ampl->Leg(i)->SetNMax(nmax);
      ampl->Leg(i)->SetKT2(0,0.);
      ampl->Leg(i)->SetKT2(1,0.);
    }
  }
  ClusterAmplitude_Vector ampls;
  p_sproc=p_proc->Caller()->Get<Single_Process>();
  ampl->SetLKF(1.0);
  int mm(p_proc->Caller()->Generator()->SetMassMode(m_nproc?0:1));
  if (!m_nproc && p_proc->Caller()->Generator()->ShiftMasses(ampl)!=1) {
    p_proc->Caller()->Generator()->SetMassMode(mm);
    ampl->Delete();
    return sqrt(-1);
  }
  Cluster(ampl,ampls,1);
  p_proc->Caller()->Generator()->SetMassMode(mm);
  if (ampls.empty()) {
    SetCoreScale(ampl);
    if (m_nproc) ampl->SetOrderQCD(ampl->OrderQCD()-1);
    p_cs->SetColors(ampl);
    if (m_nproc) ampl->SetOrderQCD(ampl->OrderQCD()+1);
    m_ampls.push_back(ampl);
    msg_Debugging()<<"Rescue: "<<*ampl<<"\n";
    return SetScales(m_ampls.back());
  }
  ampl->Delete();
  int maxlen(-1);
  std::vector<int> len(ampls.size(),-1);
  for (size_t i(0);i<ampls.size();++i) {
    for (Cluster_Amplitude *campl(ampls[i]);
	 campl;campl=campl->Next()) ++len[i];
    if (len[i]>maxlen) maxlen=len[i];
  }
  size_t imax(0);
  double sum(0.0), max(0.0);
  msg_Debugging()<<"Final configurations: max "<<maxlen<<" {\n";
  for (size_t i(0);i<ampls.size();++i) {
    msg_Indent();
    for (Cluster_Amplitude *campl(ampls[i]);
	 campl;campl=campl->Next())
      msg_Debugging()<<i<<": "<<*campl<<"\n";
    if (len[i]<maxlen) continue;
    sum+=dabs(ampls[i]->Last()->LKF());
    if (dabs(ampls[i]->Last()->LKF())>max) {
      max=dabs(ampls[i]->Last()->LKF());
      imax=i;
    }
  }
  msg_Debugging()<<"}\n";
  bool usemax((m_cmode&4) || ((m_cmode&32) && m_nproc));
  double disc(sum*ran->Get());
  sum=0.0;
  for (size_t i(0);i<ampls.size();++i) {
    if (len[i]<maxlen) continue;
    if (usemax?i==imax:(sum+=dabs(ampls[i]->Last()->LKF()))>=disc) {
      m_ampls.push_back(ampls[i]);
      ampls[i]=NULL;
      if (m_nproc) m_ampls.back()->SetOrderQCD
		     (m_ampls.back()->OrderQCD()-1);
      ClusterAmplitude_Vector ampls(m_ampls);
      m_ampls.clear();
      p_cs->SetColors(ampls.back()->Last());
      m_ampls=ampls;
      if (m_nproc) m_ampls.back()->SetOrderQCD
		     (m_ampls.back()->OrderQCD()+1);
      msg_Debugging()<<"Selected configuration "<<i<<": {\n";
      msg_Indent();
      for (Cluster_Amplitude *campl(m_ampls.back());
	   campl;campl=campl->Next()) {
	campl->SetLKF(1.0);
	msg_Debugging()<<*campl<<"\n";
      }
      break;
    }
  }
  msg_Debugging()<<"}\n";
  if (m_ampls.empty()) {
    msg_Error()<<METHOD<<"("<<p_proc->Name()<<"): Invalid point.\n";
    return -1.0;
  }
  for (size_t i(0);i<ampls.size();++i)
    if (ampls[i]) ampls[i]->Delete();
  return SetScales(m_ampls.back());
}

bool MEPS_Scale_Setter::CoreCandidate(Cluster_Amplitude *const ampl) const
{
  return ampl->Legs().size()==ampl->NIn()+m_nmin ||
    (ampl->Legs().size()==ampl->NIn()+2 &&
     ampl->Leg(2)->Flav().Mass()==0.0 &&
     ampl->Leg(3)->Flav().Mass()==0.0);
}

void MEPS_Scale_Setter::Cluster
(Cluster_Amplitude *ampl,ClusterAmplitude_Vector &ampls,const int ord) const
{
  ampl->SetProc(p_proc);
  ampl->SetProcs(p_proc->AllProcs());
  ampl->SetMS(p_proc->Generator());
  size_t oldsize(ampls.size());
  bool frs(m_rproc && ampl->Prev()==NULL);
  bool strict(!((m_cmode&1) && !rpa->gen.NumberOfTrials()));
  DEBUG_FUNC("nmin = "<<m_nmin<<", strict = "<<strict);
  msg_Debugging()<<*ampl<<"\n";
  ClusterInfo_Vector ccs;
  if (!CoreCandidate(ampl)) {
    for (size_t i(0);i<ampl->Legs().size();++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(0);j<ampl->Legs().size();++j) {
	if (i==j || (i<p_proc->NIn()&&j<p_proc->NIn())) continue;
	Cluster_Leg *lj(ampl->Leg(j));
	if (!p_sproc->Combinable(li->Id(),lj->Id())) continue;
	const Flavour_Vector &cf(p_sproc->CombinedFlavour(li->Id()+lj->Id()));
	for (size_t f(0);f<cf.size();++f) {
	  for (size_t k(0);k<ampl->Legs().size();++k) {
	    Cluster_Leg *lk(ampl->Leg(k));
	    if (k!=i && k!=j) {
	      Cluster_Config cc(ampl,i,j,k,cf[f],p_ms,NULL,-1,m_nproc?16:0);
	      if (frs && !CheckSubEvents(cc)) continue;
	      DEBUG_FUNC("Combine "<<ID(li->Id())<<" & "<<ID(lj->Id())
			 <<" <-> "<<ID(lk->Id())<<" ["<<cc.m_mo<<"], f = "<<f);
	      if (!ampl->CheckColors(li,lj,lk,cc.m_mo)) {
		msg_Debugging()<<"Veto colors: "<<li->Col()<<" & "
			       <<lj->Col()<<" <-> "<<lk->Col()<<"\n";
		continue;
	      }
	      Cluster_Info ci(cc,strict?p_clu->Cluster(cc):
			      (cc.PureQCD()?p_qdc->Cluster(cc):NULL));
	      if (ci.second.m_kt2<0.0) continue;
	      ccs.push_back(ci);
	      if (ci.second.p_ca==NULL) break;
	    }
	  }
	}
      }
    }
  }
  msg_Debugging()<<"Combinations: {\n";
  {
    msg_Indent();
    msg_Debugging()<<*ampl<<"\n";
    for (size_t i(0);i<ccs.size();++i)
      if (ccs[i].second.m_op)
	msg_Debugging()<<ccs[i].first<<" "<<ccs[i].second<<"\n";
  }
  msg_Debugging()<<"}\n";
  Int_Vector on(ccs.size(),1);
  if (frs && (m_cmode&64)) {
    msg_Debugging()<<"First RS step special\n";
    int imax(Select(ccs,on,1));
    if (imax<0) return;
    ccs[imax].second.m_op=1.0;
    if (ClusterStep(ampl,ampls,ccs[imax],0)) return;
  }
  else if (m_cmode&2) {
    for (size_t i(0);i<ccs.size();++i)
      if (ccs[i].second.m_op) ClusterStep(ampl,ampls,ccs[i],ord);
  }
  else {
    for (int i(Select(ccs,on));i>=0;on[i]=0,i=Select(ccs,on))
      if (ClusterStep(ampl,ampls,ccs[i],ord)) return;
  }
  if ((m_cmode&16) && !CoreCandidate(ampl)) return;
  if (ampls.size()>oldsize) return;
  SetCoreScale(ampl);
  if (!CheckOrdering(ampl,ord)) return;
  ampl->SetLKF((ampl->Prev()?ampl->Prev()->LKF():1.0)*Differential(ampl,1));
  if (ampl->LKF()) ampls.push_back(ampl->First()->CopyAll());
}

bool MEPS_Scale_Setter::ClusterStep
(Cluster_Amplitude *ampl,ClusterAmplitude_Vector &ampls,
 const Cluster_Info &ci,const int ord) const
{
  ampl->SetKT2(ci.second.m_kt2);
  ampl->SetMu2(ci.second.m_mu2>0.0?ci.second.m_mu2:ci.second.m_kt2);
  if (!CheckSplitting(ci,ord)) return false;
  ampl->SetLKF((ampl->Prev()?ampl->Prev()->LKF():1.0)*ci.second.m_op);
  ampl=ampl->InitNext();
  ampl->CopyFrom(ampl->Prev());
  ampl->SetCA(ci.second.p_ca);
  Combine(*ampl,ci);
  size_t oldsize(ampls.size());
  Cluster(ampl,ampls,ord);
  ampl=ampl->Prev();
  ampl->DeleteNext();
  return ampls.size()>oldsize;
}

double MEPS_Scale_Setter::Differential
(Cluster_Amplitude *const ampl,const int mode) const
{
  if (m_cmode&1024) return 1.0;
  if (ampl->Prev()==NULL) return 1.0;
  NLOTypeStringProcessMap_Map *procs
    (ampl->Procs<NLOTypeStringProcessMap_Map>());
  if (procs==NULL) return 1.0;
  nlo_type::code type=nlo_type::lo;
  if (procs->find(type)==procs->end()) {
    if (m_rproc && ampl->Prev() &&
	ampl->Prev()->Prev()==NULL) return 1.0;
    return 0.0;
  }
  Cluster_Amplitude *campl(ampl->Copy());
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  Process_Base::SortFlavours(campl);
  std::string pname(Process_Base::GenerateName(campl));
  StringProcess_Map::const_iterator pit((*(*procs)[type]).find(pname));
  if (pit==(*(*procs)[type]).end()) {
    (*(*procs)[type])[pname]=NULL;
    pit=(*procs)[type]->find(pname);
  }
  if (pit->second==NULL) {
    campl->Delete();
    return 0.0;
  }
  int kfon(pit->second->KFactorSetter(true)->On());
  pit->second->KFactorSetter(true)->SetOn(false);
  double meps = static_cast<double>(pit->second->Differential(
      *campl, Variations_Mode::nominal_only, 2 | 4 | 128 | mode));
  pit->second->KFactorSetter(true)->SetOn(kfon);
  msg_Debugging()<<"ME = "<<meps<<"\n";
  campl->Delete();
  return meps;
}

void MEPS_Scale_Setter::SetCoreScale(Cluster_Amplitude *const ampl) const
{
  PDF::Cluster_Param kt2(p_core->Calculate(ampl));
  ampl->SetKT2(kt2.m_kt2);
  ampl->SetMu2(kt2.m_mu2);
  for (Cluster_Amplitude *campl(ampl);
       campl;campl=campl->Prev()) campl->SetMuQ2(kt2.m_op);
}

double MEPS_Scale_Setter::UnorderedScale(Cluster_Amplitude *const ampl) const
{
  if (p_uoscale==NULL) return 0.;
  ampl->SetProc(p_proc);
  PDF::Cluster_Param kt2(p_uoscale->Calculate(ampl));
  return kt2.m_kt2;
}

double MEPS_Scale_Setter::SetScales(Cluster_Amplitude *ampl)
{
  double muf2(ampl->Last()->KT2()), mur2(m_rsf*ampl->Last()->Mu2());
  m_scale[stp::size+stp::res]=m_scale[stp::res]=ampl->MuQ2();
  if (ampl) {
    m_scale[stp::size+stp::res]=ampl->KT2();
    double muuo2(UnorderedScale(ampl));
    std::vector<double> scale(p_proc->NOut()+1);
    msg_Debugging()<<"Setting scales {\n";
    mur2=1.0;
    double as(1.0), mmur2(1.0), mas(1.0), oqcd(0.0), mup2(1.0);
    for (size_t idx(2);ampl->Next();++idx,ampl=ampl->Next()) {
      scale[idx]=Max(ampl->Mu2(),m_rsf*MODEL::as->CutQ2());
      scale[idx]=Min(scale[idx],sqr(rpa->gen.Ecms()));
      bool skip(false);
      Cluster_Amplitude *next(ampl->Next());
      if (!skip && next->Decays().size()) {
	size_t cid(0);
	for (size_t i(0);i<next->Legs().size();++i)
	  if (next->Leg(i)->K()) {
	    cid=next->Leg(i)->Id();
	    break;
	  }
	for (size_t i(0);i<next->Decays().size();++i)
	  if ((next->Decays()[i]->m_id&cid)==cid) {
	    skip=true;
	    break;
	  }
      }
      if (skip) continue;
      if (!s_allowuo) mup2=Max(mup2,scale[idx]);
      if (m_rproc && ampl->Prev()==NULL) {
	m_scale[stp::size+stp::res]=ampl->Next()->KT2();
	ampl->SetNLO(1);
	continue;
      }
      double coqcd(ampl->OrderQCD()-ampl->Next()->OrderQCD());
      if (coqcd>0.0) {
	if (idx>0 && scale[idx]<mup2) {
	  msg_Debugging()<<"  unordered history, redefine k_{T,"
			 <<idx<<"} = "<<sqrt(scale[idx])
			 <<" -> k_{T,"<<idx<<"} = "
			 <<sqrt(p_uoscale?muuo2:mup2)<<"\n";
	  if (p_uoscale) scale[idx]=muuo2;
	  else scale[idx]=mup2;
	}
	double cas(MODEL::as->BoundedAlphaS(m_rsf*scale[idx]));
	msg_Debugging()<<"  \\mu_{"<<idx<<"} = "
		       <<sqrt(m_rsf)<<" * "<<sqrt(scale[idx])
		       <<", as = "<<cas<<", O(QCD) = "<<coqcd<<"\n";
	mur2*=pow(m_rsf*scale[idx],coqcd);
	as*=pow(cas,coqcd);
	mmur2=Max(mmur2,m_rsf*scale[idx]);
	mas=Min(mas,cas);
	oqcd+=coqcd;
      }
      else {
        msg_Debugging()<<"  \\mu_{"<<idx<<"} = "
                       <<sqrt(m_rsf)<<" * "<<sqrt(scale[idx])
                       <<", EW splitting\n";
      }
      if (oqcd==0) m_scale[stp::size+stp::res]=ampl->Next()->KT2();
    }
    m_scale[stp::res]=ampl->MuQ2();
    double mu2(Max(ampl->Mu2(),m_rsf*MODEL::as->CutQ2()));
    double cas(MODEL::as->BoundedAlphaS(m_rsf*mu2));
    mmur2=Max(mmur2,m_rsf*mu2);
    mas=Min(mas,cas);
    if (ampl->OrderQCD()-(m_vproc?1:0)) {
      if (mu2<mup2) {
	msg_Debugging()<<"  unordered history, redefine k_{T,0} = "
		       <<sqrt(mu2)<<" -> k_{T,0} = "
		       <<sqrt(p_uoscale?muuo2:mup2)<<"\n";
	if (p_uoscale) mu2=muuo2;
	else mu2=mup2;
      }
      int coqcd(ampl->OrderQCD()-(m_vproc?1:0));
      msg_Debugging()<<"  \\mu_{0} = "<<sqrt(m_rsf)<<" * "<<sqrt(mu2)
		     <<", as = "<<cas<<", O(QCD) = "<<coqcd<<"\n";
      mur2*=pow(m_rsf*mu2,coqcd);
      as*=pow(cas,coqcd);
      oqcd+=coqcd;
    }
    if (oqcd==0) mur2=m_rsf*ampl->Mu2();
    else {
      mur2=pow(mur2,1.0/oqcd);
      as=pow(as,1.0/oqcd);
      if (m_nproc) {
	if (s_nlocpl&2) {
	  msg_Debugging()<<"  as_{NLO} = "<<mas<<"\n";
	  mur2=pow(pow(mur2,oqcd)*mmur2,1.0/(oqcd+1.0));
	  as=pow(pow(as,oqcd)*mas,1.0/(oqcd+1.0));
	}
      }
      if (s_nlocpl&1) {
	double smur2(mur2);
	mur2=MODEL::as->WDBSolve(as,m_rsf*MODEL::as->CutQ2(),
				 m_rsf*1.01*sqr(rpa->gen.Ecms()));
	if (!IsEqual(smur2,mur2))
	  msg_Debugging()<<"\\mu_R = "<<sqrt(smur2)<<" -> "<<sqrt(mur2)
			 <<", rel. dev. "<<2.*(smur2-mur2)/(smur2+mur2)<<"\n";
	if (!IsEqual((*MODEL::as)(mur2),as))
	  msg_Error()<<METHOD<<"(): Failed to determine \\mu."<<std::endl;
      }
    }
    msg_Debugging()<<"} -> as = "<<as<<" -> "<<sqrt(mur2)<<"\n";
  }
  m_scale[stp::size+stp::fac]=m_scale[stp::fac]=m_fsf*muf2;
  m_scale[stp::size+stp::ren]=m_scale[stp::ren]=mur2;
  msg_Debugging()<<"Core / QCD scale = "<<sqrt(m_scale[stp::fac])
		 <<" / "<<sqrt(m_scale[stp::ren])<<"\n";
  for (size_t i(0);i<m_calcs.size();++i)
    m_scale[i]=m_calcs[i]->Calculate()->Get<double>();
  for (size_t i(m_calcs.size());i<stp::size;++i) m_scale[i]=m_scale[0];
  if (ampl==NULL || ampl->Prev()==NULL)
    m_scale[stp::size+stp::res]=m_scale[stp::res];
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n"
		 <<"  \\mu_q = "<<sqrt(m_scale[stp::res])<<"\n";
  for (size_t i(stp::size);i<m_scale.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<(p_proc->Caller()?p_proc->Caller()->Name():"")<<"\n";
  if (ampl) {
    ampl->SetMuF2(m_scale[stp::fac]);
    ampl->SetMuR2(m_scale[stp::ren]);
    ampl->SetMuQ2(m_scale[stp::res]);
    for (size_t i(0);i<ampl->Legs().size();++i) {
      ampl->Leg(i)->SetKT2(0,-1.);
      ampl->Leg(i)->SetKT2(1,-1.);
    }
    while (ampl->Prev()) {
      ampl=ampl->Prev();
      ampl->SetMuF2(m_scale[stp::fac]);
      ampl->SetMuR2(m_scale[stp::ren]);
      ampl->SetMuQ2(m_scale[stp::res]);
      for (size_t i(0);i<ampl->Legs().size();++i) {
	ampl->Leg(i)->SetKT2(0,-1.);
	ampl->Leg(i)->SetKT2(1,-1.);
      }
    }
  }
  return m_scale[stp::fac];
}

void MEPS_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Caller()->Name()<<"' {\n";
  msg_Indent();
  m_tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}
