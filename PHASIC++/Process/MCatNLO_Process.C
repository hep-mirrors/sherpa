#include "PHASIC++/Process/MCatNLO_Process.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace PHASIC;
using namespace PDF;

MCatNLO_Process::MCatNLO_Process
(ME_Generators& gens,NLOTypeStringProcessMap_Map *pmap):
  m_gens(gens), p_bviproc(NULL), p_rsproc(NULL),
  p_bproc(NULL), p_rproc(NULL),
  p_nlomc(NULL), p_ampl(NULL)
{
  m_tinfo=0;
  static bool ref(false);
  p_apmap=pmap;
  if (!ref) {
    ref=true;
    rpa->gen.AddCitation
      (1,"The automation of MCatNLO is published in \\cite{Hoeche:2011fd}.");
  }
  m_wassevent=false;
}

MCatNLO_Process::~MCatNLO_Process()
{
  if (p_ampl) p_ampl->Delete();
  if (p_rproc) delete p_rproc;
  if (p_bproc) delete p_bproc;
  if (p_rsproc) delete p_rsproc;
  if (p_bviproc) delete p_bviproc;
}

void MCatNLO_Process::Init(const Process_Info &pi,
			  BEAM::Beam_Spectra_Handler *const beam,
			  PDF::ISR_Handler *const isr)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlo_type::born|nlo_type::loop|
		      nlo_type::vsub|nlo_type::real|nlo_type::rsub);
  if (cpi.m_fi.m_nloewtype&nlo_type::real)
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
  else
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
  Process_Base::Init(cpi,beam,isr);
  m_pinfo.m_fi.SetNLOType(pi.m_fi.NLOType());
  Process_Info npi(m_pinfo);
  npi.m_fi.m_ps.pop_back();
  m_name=GenerateName(m_pinfo.m_ii,npi.m_fi);
  if (pi.Has(nlo_type::real)!=pi.Has(nlo_type::rsub))
    THROW(fatal_error, "R/S can't be initialised separately.");
  Process_Info spi(pi);
  ++spi.m_fi.m_nmax;
  spi.m_fi.SetNLOType(cpi.m_fi.NLOType());
  p_bproc=InitProcess(spi,nlo_type::lo,false);
  p_rproc=InitProcess(spi,nlo_type::lo,true);
  p_bviproc=InitProcess(spi,nlo_type::born|nlo_type::loop|nlo_type::vsub,false);
  p_rsproc=InitProcess(spi,nlo_type::real|nlo_type::rsub,true);
  p_bviproc->FillProcessMap(p_apmap);
  p_rsproc->FillProcessMap(p_apmap);
  p_bviproc->SetSProc(p_rsproc);
  p_bproc->SetLookUp(false);
  p_rproc->SetLookUp(false);
  p_bproc->SetParent(this);
  p_rproc->SetParent(this);
  p_bproc->FillProcessMap(p_apmap);
  p_rproc->FillProcessMap(p_apmap);
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(rpa->GetPath());
  read.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
  if (!read.ReadFromFile(m_fomode,"PP_FOMODE")) m_fomode=0;
  else msg_Info()<<METHOD<<"(): Set fixed order mode "<<m_fomode<<".\n";
  if (!m_fomode) {
    p_bviproc->SetMCMode(1);
    p_rsproc->SetMCMode(2);
  }
  if (p_rsproc->Size()!=p_rproc->Size())
    THROW(fatal_error,"R and RS have different size");
  if (p_bproc->Size()!=p_bviproc->Size())
    THROW(fatal_error,"B and BVI have different size");
  for (size_t i(0);i<p_rsproc->Size();++i)
    if ((*p_rsproc)[i]->Flavours()!=(*p_rproc)[i]->Flavours())
      THROW(fatal_error,"Ordering differs in R and RS");
  for (size_t i(0);i<p_bviproc->Size();++i)
    if ((*p_bviproc)[i]->Flavours()!=(*p_bproc)[i]->Flavours())
      THROW(fatal_error,"Ordering differs in B and BVI");
}

Process_Base* MCatNLO_Process::InitProcess
(const Process_Info &pi,nlo_type::code nlotype,const bool real)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlotype);
  if (real) {
    if (cpi.m_fi.m_nloqcdtype==nlotype)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
    else if (cpi.m_fi.m_nloewtype==nlotype)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
    else THROW(fatal_error, "Internal error.");
  }
  Process_Base *proc(m_gens.InitializeProcess(cpi,false));
  if (proc==NULL) {
    msg_Error()<<cpi<<std::endl;
    THROW(not_implemented, "Process not found.");
  }
  return proc;
}

bool MCatNLO_Process::InitSubtermInfo()
{
  for (size_t i(0);i<p_rsproc->Size();++i) {
    NLO_subevtlist *subs((*p_rsproc)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      for (size_t ij(0);ij<sub->m_n;++ij)
	for (size_t k(0);k<sub->m_n;++k)
	  if (k!=ij && sub->p_fl[k]==sub->p_fl[sub->m_kt] && 
	      sub->p_fl[ij]==sub->p_fl[sub->m_ijt]) {
	    m_iinfo[sub->m_pname].insert(IDip_ID(ij,k));
	  }
      m_dinfo[subs->back()->m_pname].insert(*sub);
      std::string name(sub->Proc<Process_Base>()->Name());
      (*(*p_apmap)[nlo_type::rsub])[name]=sub->Proc<Process_Base>();
    }
  }
  return true;
}

Process_Base *MCatNLO_Process::FindProcess
(const Cluster_Amplitude *ampl,const nlo_type::code type,
 const bool error) const
{
  std::string name(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit(p_apmap->find(type)->second->find(name));
  if (pit!=p_apmap->find(type)->second->end()) return pit->second;
  if (error)
    THROW(fatal_error,"Process '"+name+"'("+ToString(type)+") not found");
  return NULL;
}

Process_Base *MCatNLO_Process::FindProcess
(const NLO_subevt *sub,const nlo_type::code type) const
{
  StringProcess_Map::const_iterator pit
    (p_apmap->find(type)->second->find(sub->m_pname));
  if (pit!=p_apmap->find(type)->second->end()) return pit->second;
  THROW(fatal_error,"Process '"+sub->m_pname+"'("+ToString(type)+") not found");
  return NULL;
}

Cluster_Amplitude *MCatNLO_Process::CreateAmplitude(const NLO_subevt *sub) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(m_nin);
  ampl->SetMS(p_rsproc->Generator());
  ampl->SetMuF2(sub->m_mu2[stp::fac]);
  ampl->SetMuR2(sub->m_mu2[stp::ren]);
  Int_Vector ci(sub->m_n,0), cj(sub->m_n,0);
  for (size_t i=0;i<sub->m_n;++i) {
    ampl->CreateLeg(i<m_nin?-sub->p_mom[i]:sub->p_mom[i],
		    i<m_nin?sub->p_fl[i].Bar():sub->p_fl[i],
		    ColorID(ci[i],cj[i]),sub->p_id[i]);
    if (!sub->IsReal() && sub->p_id[i]&(1<<sub->m_i)) {
      if ((sub->p_id[i]&(1<<sub->m_j))==0)
	THROW(fatal_error,"Internal error");
      ampl->Legs().back()->SetK(1<<sub->m_k);
    }
  }
  ampl->Decays()=*sub->p_dec;
  return ampl;
}

double MCatNLO_Process::Differential(const Vec4D_Vector &p)
{
  THROW(fatal_error,"Invalid function call");
  return m_last;
}

double MCatNLO_Process::LocalKFactor(const Cluster_Amplitude &ampl)
{
  DEBUG_FUNC(Name());
  Cluster_Amplitude *rampl(ampl.Prev());
  if (rampl->Legs().size()!=m_nin+m_nout) return 0.0;
  msg_Debugging()<<*rampl<<"\n";
  int rm(rampl->Leg(0)->Mom()[3]>0.0?1024:0);
  Process_Base *rsproc(FindProcess(rampl,nlo_type::rsub,false));
  if (rsproc==NULL) return 0.0;
  msg_Debugging()<<"Found '"<<rsproc->Name()<<"'\n";
  double rs(rsproc->Differential(*rampl,2|4|rm));
  double r(rsproc->GetSubevtList()->back()->m_result);
  msg_Debugging()<<"H = "<<rs<<", R = "<<r<<" -> "<<rs/r<<"\n";
  if (IsEqual(rs,r,1.0e-6) || r==0.0) return 1.0;
  msg_Debugging()<<ampl<<"\n";
  rm=ampl.Leg(0)->Mom()[3]>0.0?1024:0;
  Process_Base *bviproc(FindProcess(&ampl,nlo_type::vsub,false));
  if (bviproc==NULL) return 0.0;
  msg_Debugging()<<"Found '"<<bviproc->Name()<<"'\n";
  Process_Base *bproc(FindProcess(&ampl));
  double b(bproc->Differential(ampl,4|rm));
  if (b==0.0) return 0.0;
  bviproc->BBarMC()->GenerateEmissionPoint(ampl,rm);
  size_t mcmode(bviproc->SetMCMode(3));
  double bvi(bviproc->Differential(ampl,4|rm));
  bviproc->SetMCMode(mcmode);
  double k(bvi/b*(1.0-rs/r)+rs/r);
  msg_Debugging()<<"BVI = "<<bvi<<", B = "<<b
		 <<" -> K = "<<k<<"\n";
  // if (dabs(k)-1.0>r/b*ampl.Prev()->KT2()/ampl.KT2()) return 1.0;
  return k;
}

Cluster_Amplitude *MCatNLO_Process::GetAmplitude()
{
  Cluster_Amplitude *ampl(p_ampl);
  p_ampl=NULL;
  return ampl;
}

double MCatNLO_Process::OneHEvent(const int wmode)
{
  msg_Debugging()<<"H event\n";
  m_wassevent=false;
  if (p_ampl) p_ampl->Delete();
  p_selected=p_rproc;
  Process_Base *rproc(NULL);
  for (size_t i(0);i<p_rsproc->Size();++i)
    if ((*p_rsproc)[i]==p_rsproc->Selected()) {
      p_rproc->SetSelected(rproc=(*p_rproc)[i]);
      break;
    }
  rproc->Integrator()->SetMax
    (p_rsproc->Selected()->Integrator()->Max());
  Vec4D_Vector &p(p_rsproc->Selected()->Integrator()->Momenta());
  rproc->SetFixedScale(p_rsproc->Selected()->ScaleSetter(1)->Scales());
  rproc->ScaleSetter(1)->CalculateScale(Vec4D_Vector());
  rproc->SetFixedScale(std::vector<double>());
  rproc->GetMEwgtinfo()->m_mur2=
    p_rsproc->Selected()->GetMEwgtinfo()->m_mur2;
  rproc->Trigger(p);
  rproc->Integrator()->SetMomenta(p);
  p_ampl = dynamic_cast<Single_Process*>
	 (p_rsproc->Selected())->Cluster(512|2048|4096);
  if (p_ampl==NULL) {
    rproc->Differential(p);
    p_ampl = dynamic_cast<Single_Process*>(rproc)->Cluster(256|512|4096);
  }
  if (p_ampl==NULL) {
    msg_Error()<<METHOD<<"(): No valid clustering. Skip event."<<std::endl;
    return 0.0;
  }
  Scale_Setter_Base *scs(p_rsproc->Selected()->ScaleSetter(1));
  for (Cluster_Amplitude *ampl(p_ampl);ampl;ampl=ampl->Next()) {
    ampl->SetMuR2(scs->Scale(stp::ren));
    ampl->SetMuF2(scs->Scale(stp::fac));
    ampl->SetQ2(scs->Scale(stp::res));
  }
  Selector_Base *jf=(*p_bviproc)[0]->
    Selector()->GetSelector("Jetfinder");
  if (jf && m_nout-1<m_pinfo.m_fi.NMaxExternal()) {
    bool res(p_shower->JetVeto(p_ampl));
    if (res) return 0.0;
  }
  return 1.0;
}

double MCatNLO_Process::OneSEvent(const int wmode)
{
  msg_Debugging()<<"S event\n";
  m_wassevent=true;
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;
  Process_Base *bproc(NULL);
  for (size_t i(0);i<p_bviproc->Size();++i)
    if ((*p_bviproc)[i]==p_bviproc->Selected()) {
      p_bproc->SetSelected(bproc=(*p_bproc)[i]);
      break;
    }
  Vec4D_Vector &p(p_bviproc->Selected()->Integrator()->Momenta());
  bproc->Trigger(p);
  bproc->Differential(p);
  p_ampl = dynamic_cast<Single_Process*>(bproc)->Cluster(256|512|4096);
  SortFlavours(p_ampl);
  p_ampl->SetProcs(p_apmap);
  p_ampl->SetIInfo(&m_iinfo);
  p_ampl->SetDInfo(&m_dinfo);
  p_ampl->Decays()=m_decins;
  p_nlomc->SetShower(p_shower);
  int stat(p_nlomc->GeneratePoint(p_ampl));
  Cluster_Amplitude *next(p_ampl), *ampl(p_ampl->Prev());
  if (ampl) {
    p_ampl=NULL;
    Process_Base *rproc(FindProcess(ampl));
    if (rproc==NULL) THROW(fatal_error,"Invalid splitting");
    p_selected=p_rproc;
    p_rproc->SetSelected(rproc);
    rproc->Integrator()->SetMax(bproc->Integrator()->Max());
    if (ampl->Leg(0)->Mom().PPlus()>ampl->Leg(1)->Mom().PPlus())
      std::swap<Cluster_Leg*>(ampl->Legs()[0],ampl->Legs()[1]);
    rproc->Integrator()->SetMomenta(*ampl);
    rproc->GetMEwgtinfo()->m_mur2=bproc->GetMEwgtinfo()->m_mur2;
    Cluster_Leg *lij(NULL);
    for (size_t i(0);i<next->Legs().size();++i)
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
	break;
      }
    std::vector<int> ids(ID(lij->Id()));
    size_t iid(1<<ids.front()), jid(1<<ids.back()), kid(lij->K());
    ids.push_back(ID(lij->K()).front());
    if (ids.size()!=3) THROW(fatal_error,"Internal error");
    Cluster_Amplitude *campl(ampl->Copy());
    campl->IdSort();
    Cluster_Definitions_Base *clus(p_shower->GetClusterDefinitions());
    CParam kt2(clus->KPerp2(*campl,ids[0],ids[1],ids[2],
			    lij->Flav(),p_bviproc->Generator(),1));
    campl->Delete();
    ampl->SetKT2(kt2.m_kt2);
    ampl->SetMu2(kt2.m_kt2);
    p_ampl=ampl;
    ampl->SetMuF2(next->MuF2());
    ampl->SetMuR2(next->MuR2());
    ampl->SetQ2(kt2.m_kt2);
    ampl->SetOrderQCD(next->OrderQCD()+1);
    ampl->Next()->SetNLO(4);
    next->SetKin(kt2.m_kin);
    while (ampl->Next()) {
      ampl=ampl->Next();
      ampl->SetQ2(kt2.m_kt2);
      if (!(ampl->Leg(0)->Id()&p_ampl->Leg(0)->Id()))
	std::swap<Cluster_Leg*>(ampl->Legs()[0],ampl->Legs()[1]);
      if (ampl->IdNew()&iid) ampl->SetIdNew(ampl->IdNew()|jid);
      for (size_t i(0);i<ampl->Legs().size();++i) {
	ampl->Leg(i)->SetNMax(p_ampl->Leg(i)->NMax());
	size_t cid(ampl->Leg(i)->Id());
	if (cid&iid) {
	  for (size_t j(0);j<ampl->Legs().size();++j)
	    if (ampl->Leg(j)->K()==cid)
	      ampl->Leg(j)->SetK(cid|jid);
	  ampl->Leg(i)->SetId(cid|jid);
	  if (ampl->Prev()->Prev()==NULL) {
	    ampl->Leg(i)->SetK(kid);
	  }
	}
      }
    }
    msg_Debugging()<<"R selected via Sudakov "<<*p_ampl
		   <<" ( w = "<<p_nlomc->Weight()<<" )\n";
    return p_nlomc->Weight();
  }
  p_selected=p_bproc;
  if (p_ampl->Leg(0)->Mom().PPlus()>p_ampl->Leg(1)->Mom().PPlus())
    std::swap<Cluster_Leg*>(p_ampl->Legs()[0],p_ampl->Legs()[1]);
  ampl=p_ampl;
  ampl->SetNLO(4);
  bproc->Integrator()->SetMomenta(*p_ampl);
  msg_Debugging()<<"B selected "<<*p_ampl
		 <<" ( w = "<<p_nlomc->Weight()<<" )\n";
  return stat?p_nlomc->Weight():0.0;
}

Weight_Info *MCatNLO_Process::OneEvent(const int wmode,const int mode)
{
  DEBUG_FUNC("");
  double S(p_bviproc->Integrator()->SelectionWeight(wmode));
  double H(p_rsproc->Integrator()->SelectionWeight(wmode));
  Weight_Info *winfo(NULL);
  if (S>ran->Get()*(S+H)) {
    p_selected=p_bviproc;
    winfo=p_bviproc->OneEvent(wmode,mode);
    if (winfo && m_fomode==0) {
      winfo->m_weight*=OneSEvent(wmode);
      winfo->m_weight*=p_bproc->Selected()
	->Integrator()->SelectionWeight(wmode)/
	p_bviproc->Selected()->Integrator()
	->SelectionWeight(wmode);
      double lkf(p_bviproc->Selected()->Last()/
		 p_bproc->Selected()->Last());
      for (Cluster_Amplitude *ampl(p_ampl);
	   ampl;ampl=ampl->Next()) ampl->SetLKF(lkf);
      Cluster_Amplitude *ampl(p_ampl->Next());
      if (ampl) {
	Selector_Base *jf=(*p_bviproc)[0]->
	  Selector()->GetSelector("Jetfinder");
	if (jf) {
	  if (ampl->NLO()!=0) ampl=ampl->Next();
	  for (;ampl;ampl=ampl->Next()) ampl->SetNLO(2);
	}
      }
    }
  }
  else {
    p_selected=p_rsproc;
    p_rsproc->SetClusterMode(1|1024);
    winfo=p_rsproc->OneEvent(wmode,mode);
    p_rsproc->SetClusterMode(0);
    if (winfo && m_fomode==0) {
      winfo->m_weight*=OneHEvent(wmode);
      winfo->m_weight*=p_rproc->Selected()
	->Integrator()->SelectionWeight(wmode)/
	p_rsproc->Selected()->Integrator()
	->SelectionWeight(wmode);
    }
  }
  for (Cluster_Amplitude *ampl(p_ampl);ampl;
       ampl=ampl->Next()) ampl->SetNLO(1|ampl->NLO());
  if (winfo && winfo->m_weight==0) {
    delete winfo;
    winfo=NULL;
  }
  return winfo;
}

void MCatNLO_Process::InitPSHandler
(const double &maxerror,const std::string eobs,const std::string efunc)
{
  p_bviproc->InitPSHandler(maxerror,eobs,efunc);
  p_rsproc->InitPSHandler(maxerror,eobs,efunc);
  p_rsproc->Integrator()->SetEnhanceFactor
    (p_rsproc->Integrator()->EnhanceFactor()*p_int->RSEnhanceFactor());
  p_rproc->Integrator()->SetPSHandler
    (p_rsproc->Integrator()->PSHandler());
  p_bproc->Integrator()->SetPSHandler
    (p_bviproc->Integrator()->PSHandler());
}

bool MCatNLO_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  Vec4D_Vector p(p_rsproc->NIn()+p_rsproc->NOut());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i(0);i<p.size();++i)
    ampl->CreateLeg(Vec4D(),Flavour(kf_jet));
  SP(Phase_Space_Handler) psh(p_rsproc->Integrator()->PSHandler());
  do {
    psh->TestPoint(&p.front(),&p_rsproc->Info(),p_rsproc->Generator());
    for (size_t i(0);i<p.size();++i) ampl->Leg(i)->SetMom(p[i]);
    p_rsproc->Differential(*ampl,4);
  } while (!InitSubtermInfo());
  ampl->Delete();
  bool res(p_bviproc->CalculateTotalXSec
	   (resultpath+"/"+p_bviproc->Generator()->Name(),create));
  if (!p_rsproc->CalculateTotalXSec
      (resultpath+"/"+p_rsproc->Generator()->Name(),create)) res=false;
  for (size_t i(0);i<p_bviproc->Size();++i)
    (*p_bproc)[i]->Integrator()->SetMax
      ((*p_bviproc)[i]->Integrator()->Max());
  return res;
}

void MCatNLO_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  p_bviproc->SetLookUp(lookup);
  p_rsproc->SetLookUp(lookup);
}

void MCatNLO_Process::SetScale(const Scale_Setter_Arguments &scale)
{
  p_bviproc->SetScale(scale);
  p_rsproc->SetScale(scale);
  p_rproc->SetScale(scale);
  p_bproc->SetScale(scale);
}

void MCatNLO_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  p_bviproc->SetKFactor(args);
  p_rsproc->SetKFactor(args);
  p_rproc->SetKFactor(args);
  p_bproc->SetKFactor(args);
}

void MCatNLO_Process::SetFixedScale(const std::vector<double> &s)
{
  p_bviproc->SetFixedScale(s);
  p_rsproc->SetFixedScale(s);
  p_rproc->SetFixedScale(s);
  p_bproc->SetFixedScale(s);
}

bool MCatNLO_Process::IsGroup() const
{
  return true;
}

size_t MCatNLO_Process::Size() const
{
  return 2;
}

Process_Base *MCatNLO_Process::operator[](const size_t &i)
{
  if (i==1) return p_rsproc;
  return p_bviproc;
}

void MCatNLO_Process::SetClusterDefinitions
(PDF::Cluster_Definitions_Base *const cluster)
{
  p_bproc->Generator()->SetClusterDefinitions(cluster);
  p_rproc->Generator()->SetClusterDefinitions(cluster);
}

void MCatNLO_Process::SetSelector(const Selector_Key &key)
{
  p_bviproc->SetSelector(key);
  p_rsproc->SetSelector(key);
  p_rproc->SetSelector(key);
  p_bproc->SetSelector(key);
}

void MCatNLO_Process::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  p_bviproc->SetShower(ps);
  p_rsproc->SetShower(ps);
  p_rproc->SetShower(ps);
  p_bproc->SetShower(ps);
}
