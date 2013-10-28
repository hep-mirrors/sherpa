#include "CSSHOWER++/Main/CS_Shower.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace CSSHOWER;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

CS_Shower::CS_Shower(PDF::ISR_Handler *const _isr,
		     MODEL::Model_Base *const model,
		     Data_Reader *const _dataread) : 
  Shower_Base("CSS"), p_isr(_isr), 
  p_shower(NULL), p_cluster(NULL), p_ampl(NULL)
{
  rpa->gen.AddCitation
    (1,"The Catani-Seymour subtraction based shower is published under \\cite{Schumann:2007mg}.");
  int maxem=_dataread->GetValue<int>("CSS_MAXEM",-1);
  if (maxem<0) m_maxem=std::numeric_limits<size_t>::max();
  else {
    m_maxem=maxem;
    msg_Info()<<METHOD<<"(): Set max emissions "<<m_maxem<<"\n";
  }
  SF_Lorentz::SetKappa(_dataread->GetValue<double>("DIPOLE_KAPPA",2.0/3.0));
  m_kmode=_dataread->GetValue<int>("CSS_KMODE",2);
  if (m_kmode!=2) msg_Info()<<METHOD<<"(): Set kernel mode "<<m_kmode<<"\n";
  m_recocheck=_dataread->GetValue<int>("CSS_RECO_CHECK",0);
  if (m_recocheck!=0) msg_Info()<<METHOD<<"(): Set reco check mode "<<m_recocheck<<"\n";
  m_respectq2=_dataread->GetValue<int>("CSS_RESPECT_Q2",1);
  if (m_respectq2!=1) msg_Info()<<METHOD<<"(): Set respect Q2 mode "<<m_respectq2<<"\n";
  int amode(_dataread->GetValue<int>("EXCLUSIVE_CLUSTER_MODE",0));
  if (amode!=0) msg_Info()<<METHOD<<"(): Set exclusive cluster mode "<<amode<<".\n";
  int meweight=_dataread->GetValue<int>("CSS_MEWMODE",1);
  if (meweight!=1) msg_Info()<<METHOD<<"(): Set ME weight mode "<<meweight<<"\n";
  
  m_weightmode = int(_dataread->GetValue<int>("WEIGHT_MODE",1));
  
  int _qed=_dataread->GetValue<int>("CSS_EW_MODE",0);
  if (_qed==1) {
    s_kftable[kf_photon]->SetResummed();
  }
  p_shower = new Shower(_isr,_qed,_dataread);
  
  p_next = new All_Singlets();

  p_cluster = new CS_Cluster_Definitions(p_shower,m_kmode,meweight);
  p_cluster->SetAMode(amode);
}

CS_Shower::~CS_Shower() 
{
  CleanUp();
  if (p_shower)      { delete p_shower; p_shower = NULL; }
  if (p_cluster)     { delete p_cluster; p_cluster = NULL; }
  if (p_ampl) p_ampl->Delete();
  delete p_next;
}

int CS_Shower::PerformShowers(const size_t &maxem,size_t &nem)
{
  if (!p_shower || !m_on) return 1;
  m_weight=1.0;
  for (All_Singlets::const_iterator sit(m_allsinglets.begin());
       sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before shower step\n";
    for (Singlet::const_iterator it((*sit)->begin());it!=(*sit)->end();++it)
      if ((*it)->GetPrev() && (*it)->GetPrev()->KScheme()!=1) {
	if (m_respectq2) (*it)->SetStart(Min((*it)->KtStart(),(*it)->GetPrev()->KtStart()));
	else (*it)->SetStart((*it)->GetPrev()->KtStart());
      }
    msg_Debugging()<<**sit;
    size_t pem(nem);
    if (!p_shower->EvolveShower(*sit,maxem,nem)) return 0;
    m_weight*=p_shower->Weight();
    if ((*sit)->GetLeft()) {
      p_shower->ReconstructDaughters(*sit,1);
    }
    msg_Debugging()<<"after shower step with "<<nem-pem
		   <<" of "<<nem<<" emission(s)\n";
    msg_Debugging()<<**sit<<"\n";
  }
  return 1;
}

int CS_Shower::PerformShowers() 
{
  return PerformShowers(m_maxem,m_nem);
}

int CS_Shower::PerformDecayShowers() {
  if (!p_shower) return 1;
  size_t nem(0);
  for (All_Singlets::const_iterator 
	 asit(m_allsinglets.begin());asit!=m_allsinglets.end();++asit) {
    if (!p_shower->EvolveShower(*asit,m_maxem,nem)) return 0;
  }
  return 1;
}

bool CS_Shower::ExtractPartons(Blob_List *const blist) {
  
  Blob * psblob(blist->FindLast(btp::Shower));
  if (psblob==NULL) THROW(fatal_error,"No Shower blob");
  psblob->SetTypeSpec("CSSHOWER++1.0");
  for (int i=0;i<psblob->NInP();++i)
    psblob->InParticle(i)->SetStatus(part_status::decayed);
  for (int i=0;i<psblob->NOutP();++i)
    psblob->OutParticle(i)->SetStatus(part_status::decayed);
  
  psblob->SetStatus(blob_status::needs_beams |
		    blob_status::needs_hadronization);
  
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit)
      (*sit)->ExtractPartons(psblob,p_ms);
  return true;
}

void CS_Shower::CleanUp()
{
  m_nem=0;
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    if (*sit) delete *sit;
  }
  m_allsinglets.clear();
}

PDF::Cluster_Definitions_Base * CS_Shower::GetClusterDefinitions() 
{
  return p_cluster;
}

void CS_Shower::GetKT2Min(Cluster_Amplitude *const ampl,const size_t &id,
			  KT2X_Map &kt2xmap,std::set<size_t> &aset)
{
  msg_Indent();
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if ((cl->Id()&id)==0) continue;
    if (ampl->Prev()) GetKT2Min(ampl->Prev(),cl->Id(),kt2xmap,aset);
    if (cl->Stat()&2) {
      double ckt2(dabs(cl->Mom().Abs2()));
      kt2xmap[cl->Id()].first=kt2xmap[cl->Id()].second=HardScale(ampl);
      for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
	if (kit->first!=cl->Id() && (kit->first&cl->Id()) &&
	    aset.find(kit->first)==aset.end()) {
	  kit->second.first=ckt2;
	  kit->second.second=ckt2;
	  aset.insert(kit->first);
	}
    }
    else if (ampl->Prev()==NULL) {
      kt2xmap[cl->Id()].first=kt2xmap[cl->Id()].second=HardScale(ampl); 
    }
    else {
      double ckt2max(HardScale(ampl));
      double ckt2min(std::numeric_limits<double>::max());
      for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
	if (kit->first&cl->Id()) {
	  ckt2min=kit->second.first;
	  ckt2max=Max(ckt2max,kit->second.second);
	}
      kt2xmap[cl->Id()].first=ckt2min;
      kt2xmap[cl->Id()].second=ckt2max; 
      for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
	if ((kit->first&cl->Id()) && aset.find(kit->first)==aset.end()) {
	  kit->second.first=ckt2min;
	  kit->second.second=ckt2max;
	}
    }
  }
}

void CS_Shower::GetKT2Min(Cluster_Amplitude *const ampl,KT2X_Map &kt2xmap)
{
  std::set<size_t> aset;
  Cluster_Amplitude *campl(ampl);
  while (campl->Next()) campl=campl->Next();
  GetKT2Min(campl,(1<<ampl->Legs().size())-1,kt2xmap,aset);
  std::vector<size_t> cns;
  double ckt2min(std::numeric_limits<double>::max()), ckt2max(0.0);
  for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
    if (aset.find(kit->first)==aset.end()) {
      ckt2min=Min(ckt2min,kit->second.first);
      ckt2max=Max(ckt2max,kit->second.second);
      bool ins(true);
      for (size_t j(0);j<cns.size();++j)
	if (cns[j]&kit->first) {
	  ins=false;
	  break;
	}
      if (ins) cns.push_back(kit->first);
    }
  Cluster_Amplitude *rampl(ampl);
  for (;rampl->Next();rampl=rampl->Next()) {
    bool dc(false);
    for (size_t i(0);i<rampl->Next()->Legs().size();++i)
      if (rampl->Next()->Leg(i)->Stat()&4) dc=true;
    if (!dc) break;
  }
  bool smin(rampl->Legs().size()-rampl->NIn()==campl->Leg(0)->NMax());
  for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
    if (aset.find(kit->first)==aset.end()) {
      if (smin) kit->second.first=ckt2min;
      else kit->second.first=0.0;
      kit->second.second=ckt2max;
    }
  msg_Debugging()<<"k_{T,min} / k_{T,max} = {\n";
  for (KT2X_Map::const_iterator
  	 kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
    msg_Debugging()<<"  "<<ID(kit->first)
		  <<" -> "<<sqrt(kit->second.first)
		  <<" / "<<sqrt(kit->second.second)<<"\n";
  msg_Debugging()<<"}\n";
}

bool CS_Shower::PrepareShower(Cluster_Amplitude *const ampl,const bool & soft)
{
  if (soft) return PrepareShowerFromSoft(ampl);
  return PrepareStandardShower(ampl);
}


bool CS_Shower::EstablishRelations(Parton * parton,Cluster_Leg * leg,
				   std::map<Cluster_Leg*,Parton*> & pmap) {
  int no = leg->NumberOfSpectators();
  if (no==0) return true;

  size_t col1(parton->GetFlow(1)),col2(parton->GetFlow(2));
  bool connect;
  for (std::list<Cluster_Leg *>::iterator clit=leg->GetSpectators().begin();
       clit!=leg->GetSpectators().end();clit++) {
    connect = false;
    Parton * spect = pmap[(*clit)];
    if (col1!=0 && col1==spect->GetFlow(2)) {
      if ((parton->GetLeft() && parton->GetLeft()!=spect) ||
	  (spect->GetRight() && spect->GetRight()!=parton)) {
	connect = true;
      }
      else {
	parton->SetLeft(spect);
	spect->SetRight(parton);
	connect = true;
      }
    }
    if (col2!=0 && col2==spect->GetFlow(1)) {
      if ((parton->GetRight() && parton->GetRight()!=spect) ||
	  (spect->GetLeft() && spect->GetLeft()!=parton)) {
	connect = true;
      }
      else {
	parton->SetRight(spect);
	spect->SetLeft(parton);
	connect = true;
      }
    }
    if (!connect) { 
      if (spect->GetLeft()==parton)       parton->SetRight(spect);
      else if (spect->GetRight()==parton) parton->SetLeft(spect);
      else {
	if (!parton->GetLeft())       parton->SetLeft(spect);
	else if (!parton->GetRight()) parton->SetRight(spect);
      }
    }
  }
  if (parton->GetFlavour().IsGluon()) {
    if (parton->GetLeft()==NULL) {
      msg_Debugging()<<METHOD<<": gluon "<<ATOOLS::ID(parton->Id())
		     <<" with no left - use "
		     <<ATOOLS::ID(parton->GetRight()->Id())<<" instead.\n";
      parton->SetLeft(parton->GetRight());
    }
    if (parton->GetRight()==NULL) {
      msg_Debugging()<<METHOD<<": gluon "<<ATOOLS::ID(parton->Id())
		     <<" with no right - use "
		     <<ATOOLS::ID(parton->GetLeft()->Id())<<" instead.\n";
      parton->SetRight(parton->GetLeft());
    }
  }
  if (parton->GetRight()) parton->GetRight()->SetInvRight(parton);
  if (parton->GetLeft())  parton->GetLeft()->SetInvLeft(parton);
  parton->SetStat(leg->Stat());
  return true;
}

bool CS_Shower::PrepareShowerFromSoft(Cluster_Amplitude *const ampl)
{
  CleanUp();
  p_rampl = ampl;
  p_ms    = ampl->MS();
  p_next->clear();
  std::map<Parton*,Cluster_Leg*,partcomp> lmap;
  std::map<Cluster_Leg*,Parton*> pmap;

  Parton      * parton;
  Cluster_Leg * leg;
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
  for (size_t i(0);i<ampl->Legs().size();++i) {
    leg = ampl->Leg(i);
    if (leg->Flav().IsHadron() && 
	leg->Flav().Kfcode()!=kf_pomeron &&
	leg->Flav().Kfcode()!=kf_reggeon && 
	leg->Id()&((1<<ampl->NIn())-1)) continue;
    bool is(leg->Id()&((1<<ampl->NIn())-1));
    Particle p(1,is?leg->Flav().Bar():leg->Flav(),is?-leg->Mom():leg->Mom());
    if (is) {
      p.SetFlow(2,leg->Col().m_i);
      p.SetFlow(1,leg->Col().m_j);
    }
    else {
      p.SetFlow(1,leg->Col().m_i);
      p.SetFlow(2,leg->Col().m_j);
    }
    parton       = new Parton(&p,is?pst::IS:pst::FS);
    pmap[leg]    = parton;
    lmap[parton] = leg;
    parton->SetRFlow();
    parton->SetKin(p_shower->KinScheme());
    parton->SetId(leg->Id());
    if (is) {
      if (Vec3D(p.Momentum())*Vec3D(rpa->gen.PBeam(0))>0.) {
	parton->SetBeam(0);
      }
      else { 
	parton->SetBeam(1);
      }
    }
    double kt2start(leg->KTStart());
    parton->SetStart(kt2start);   // start scale of shower
    parton->SetKtMax(sqr(rpa->gen.Ecms()));    // no jet veto below ktmax
    parton->SetVeto(sqr(rpa->gen.Ecms()));     // irrelevant 
    parton->SetConnected(leg->Connected());
    parton->SetMass2(p_ms->Mass2(leg->Flav()));
    singlet->push_back(parton);
    parton->SetSing(singlet);
  }
  for (std::map<Parton*,Cluster_Leg*>::iterator pit=lmap.begin();
       pit!=lmap.end();pit++) {
    parton = pit->first;
    leg    = pit->second;
    if (!EstablishRelations(parton,leg,pmap)) {
      for (size_t i(0);i<ampl->Legs().size();++i) {
	msg_Debugging()<<(*ampl->Leg(i))<<"\n";
      }
      return false;
    }
  }
  msg_Debugging()<<METHOD<<" for "<<lmap.size()<<" partons:\n";
  for (std::map<Parton*,Cluster_Leg*>::iterator pit=lmap.begin();
       pit!=lmap.end();pit++) {
    parton = pit->first;
    if (parton->GetType()==pst::IS) continue;
    msg_Debugging()<<parton<<" ("<<parton->GetFlavour()<<") "
		   <<"["<<parton->GetFlow(1)<<", "<<parton->GetFlow(2)<<"], "
		   <<"{"<<parton->GetLeft()<<", "<<parton->GetRight()<<"}; "
		   <<"{"<<parton->GetInvLeft()<<", "<<parton->GetInvRight()<<"};\n";
  }
  p_shower->SetMS(p_ms);

  pmap.clear();
  lmap.clear();
  m_allsinglets.push_back(singlet);
  return true;
}

bool CS_Shower::PrepareStandardShower(Cluster_Amplitude *const ampl)
{
  CleanUp();
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  p_rampl=ampl;
  p_ms=ampl->MS();
  KT2X_Map kt2xmap;
  GetKT2Min(ampl,kt2xmap);
  p_next->clear();
  All_Singlets allsinglets;
  std::map<size_t,Parton*> kmap;
  std::map<Parton*,Cluster_Leg*> almap;
  std::map<Cluster_Leg*,Parton*> apmap;
  for (Cluster_Amplitude *campl(ampl);campl;campl=campl->Next()) {
    msg_Debugging()<<*campl<<"\n";
    Parton *split(NULL);
    std::map<Parton*,Cluster_Leg*> lmap;
    std::map<Cluster_Leg*,Parton*> pmap;
    Singlet *sing(TranslateAmplitude(campl,pmap,lmap,kt2xmap));
    allsinglets.push_back(sing);
    for (size_t i(0);i<campl->Legs().size();++i) {
      Cluster_Leg *cl(campl->Leg(i));
      if (pmap.find(cl)==pmap.end()) continue;
      Parton *k(pmap[cl]);
      k->SetId(cl->Id());
      almap[apmap[cl]=k]=cl;
      std::map<size_t,Parton*>::iterator cit(kmap.find(cl->Id()));
      if (cit!=kmap.end()) {
	if (k->GetNext()!=NULL) 
	  THROW(fatal_error,"Invalid tree structure");
	k->SetNext(cit->second);
	cit->second->SetPrev(k);
	k->SetKScheme(cit->second->KScheme());
	if (k->KScheme()) k->SetMass2(k->Momentum().Abs2());
	k->SetStat(1);
	kmap[cl->Id()]=k;
	continue;
      }
      for (std::map<size_t,Parton*>::iterator 
	     kit(kmap.begin());kit!=kmap.end();)
	if ((kit->first&cl->Id())!=kit->first) {
	  ++kit;
	}
	else {
	  if (k->GetSing()->GetLeft()==NULL) k->GetSing()->SetLeft(kit->second);
	  else if (k->GetSing()->GetRight()==NULL) {
	    if (kit->second->GetType()==pst::IS &&
		k->GetSing()->GetLeft()->GetType()==pst::FS) {
	      k->GetSing()->SetRight(k->GetSing()->GetLeft());
	      k->GetSing()->SetLeft(kit->second);
	    }
	    else {
	      k->GetSing()->SetRight(kit->second);
	    }
	  }
	  else THROW(fatal_error,"Invalid tree structure");
	  kit->second->SetPrev(k);
	  k->SetStat(1);
	  kmap.erase(kit);
	  kit=kmap.begin();
	  split=k;
	}
      kmap[cl->Id()]=k;
    }
    if (sing->GetSpec()) {
      split->SetOldMomentum(split->Momentum());
      sing->SetSpec(sing->GetSpec()->GetNext());
      if (split==NULL) THROW(fatal_error,"Invalid tree structure");
      sing->SetSplit(split);
      Parton *l(sing->GetLeft()), *r(sing->GetRight()), *s(sing->GetSpec());
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      split->SetKin(campl->Kin());
      split->SetKScheme((almap[split]->Stat()&2)?1:0);
      if (split->KScheme()) split->SetMass2(split->Momentum().Abs2());
      CS_Parameters cp(p_cluster->KT2
		       (campl->Prev(),almap[l],almap[r],almap[s],
			split->GetType()==pst::FS?split->GetFlavour():
			split->GetFlavour().Bar(),p_ms,
			split->Kin(),split->KScheme()));
      l->SetTest(cp.m_kt2,cp.m_z,cp.m_y,cp.m_phi);
      if (split->KScheme()) split->SetFixSpec(cp.m_pk);
      msg_Debugging()<<"Set reco params: kt = "<<sqrt(cp.m_kt2)<<", z = "
		     <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi
		     <<", mode = "<<cp.m_mode<<", scheme = "<<split->Kin()
		     <<", kmode = "<<split->KScheme()<<"\n";
      sing->SetAll(p_next);
      Vec4D oldl(l->Momentum()), oldr(r->Momentum()), olds(s->Momentum());
      if (m_recocheck&1) {
      std::cout.precision(12);
      Vec4D oldfl(l->FixSpec()), oldfr(r->FixSpec()), oldfs(s->FixSpec());
      Vec4D oldsf(split->FixSpec()), oldso(split->OldMomentum());
      sing->BoostBackAllFS(l,r,s,split,split->GetFlavour(),cp.m_mode|4|8);
      p_shower->ReconstructDaughters(sing,1);
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      CS_Parameters ncp(p_cluster->KT2
			(campl->Prev(),almap[l],almap[r],almap[s],
			 split->GetType()==pst::FS?split->GetFlavour():
			 split->GetFlavour().Bar(),p_ms,
			 split->Kin(),split->KScheme()));
      msg_Debugging()<<"New reco params: kt = "<<sqrt(ncp.m_kt2)<<", z = "
		     <<ncp.m_z<<", y = "<<ncp.m_y<<", phi = "<<ncp.m_phi
		     <<", kin = "<<ncp.m_kin<<"\n";
      msg_Debugging()<<"            vs.: kt = "<<sqrt(cp.m_kt2)<<", z = "
		     <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi
		     <<", kin = "<<cp.m_kin<<"\n";
      if (!IsEqual(ncp.m_kt2,cp.m_kt2,1.0e-6) || 
	  !IsEqual(ncp.m_z,cp.m_z,1.0e-6) || 
	  !IsEqual(ncp.m_y,cp.m_y,1.0e-6) || 
	  !IsEqual(ncp.m_phi,cp.m_phi,1.0e-6) ||
	  !IsEqual(oldl,l->Momentum(),1.0e-6) || 
	  !IsEqual(oldr,r->Momentum(),1.0e-6) || 
	  !IsEqual(olds,s->Momentum(),1.0e-6)) {
	msg_Error()<<"\nFaulty reco params: kt = "<<sqrt(ncp.m_kt2)<<", z = "
		   <<ncp.m_z<<", y = "<<ncp.m_y<<", phi = "<<ncp.m_phi
		   <<", mode = "<<ncp.m_mode<<", scheme = "<<ncp.m_kin<<"\n";
	msg_Error()<<"               vs.: kt = "<<sqrt(cp.m_kt2)<<", z = "
		   <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi
		   <<", mode = "<<cp.m_mode<<", scheme = "<<cp.m_kin<<"\n";
	msg_Error()<<"  "<<oldl<<" "<<oldr<<" "<<olds<<"\n";
	msg_Error()<<"  "<<l->Momentum()<<" "<<r->Momentum()
		   <<" "<<s->Momentum()<<"\n";
	if (m_recocheck&2) abort();
      }
      l->SetMomentum(oldl);
      r->SetMomentum(oldr);
      s->SetMomentum(olds);
      l->SetFixSpec(oldfl);
      r->SetFixSpec(oldfr);
      s->SetFixSpec(oldfs);
      split->SetFixSpec(oldsf);
      split->SetOldMomentum(oldso);
      }
      l->SetOldMomentum(oldl);
      r->SetOldMomentum(oldr);
      s->SetOldMomentum(olds);
      sing->BoostBackAllFS(l,r,s,split,split->GetFlavour(),cp.m_mode|4|8);
    }
    double kt2next(campl->Prev()?campl->Prev()->KT2():0.0);
    sing->SetKtNext(kt2next);
    if (ampl->NIn()==1 && ampl->Leg(0)->Flav().IsHadron()) break;
    p_next->push_back(sing);
  }
  p_next->clear();
  for (All_Singlets::reverse_iterator
	 asit(allsinglets.rbegin());asit!=allsinglets.rend();++asit) {
    m_allsinglets.push_back(*asit);
    p_next->push_back(*asit);
  }
  msg_Debugging()<<"\nSinglet lists:\n\n";
  Cluster_Amplitude *campl(ampl);
  while (campl->Next()) campl=campl->Next();
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
      for (Singlet::const_iterator 
	     pit((*sit)->begin());pit!=(*sit)->end();++pit) {
	if ((*pit)->GetPrev()) {
	  if ((*pit)->GetPrev()->GetNext()==*pit) 
	    (*pit)->SetStart((*pit)->GetPrev()->KtStart());
	}
      }
      (*sit)->SetDecays(campl->Decays());
      (*sit)->SetAll(p_next);
      msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  msg_Debugging()<<"}\n";
  p_shower->SetMS(p_ms);
  return true;
}

Singlet *CS_Shower::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap,
 const KT2X_Map &kt2xmap)
{
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  double ktveto2(jf?jf->Ycut()*sqr(rpa->gen.Ecms()):
		 sqrt(std::numeric_limits<double>::max()));
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
  singlet->SetJF(jf);
  singlet->SetNLO(ampl->NLO()&~1);
  if (jf==NULL && (ampl->NLO()&2)) singlet->SetNLO(4);
  singlet->SetLKF(ampl->LKF());
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if (cl->Flav().IsHadron() && cl->Id()&((1<<ampl->NIn())-1)) continue;
    bool is(cl->Id()&((1<<ampl->NIn())-1));
    Particle p(1,is?cl->Flav().Bar():cl->Flav(),is?-cl->Mom():cl->Mom());
    if (is) {
      p.SetFlow(2,cl->Col().m_i);
      p.SetFlow(1,cl->Col().m_j);
    }
    else {
      p.SetFlow(1,cl->Col().m_i);
      p.SetFlow(2,cl->Col().m_j);
    }
    Parton *parton(new Parton(&p,is?pst::IS:pst::FS));
    parton->SetMass2(p_ms->Mass2(p.Flav()));
    pmap[cl]=parton;
    lmap[parton]=cl;
    parton->SetRFlow();
    parton->SetKin(p_shower->KinScheme());
    if (is) parton->SetBeam(i);
    KT2X_Map::const_iterator xit(kt2xmap.find(cl->Id()));
    parton->SetStart(m_respectq2?ampl->Q2():xit->second.second);
    parton->SetKtMax(xit->second.first);
    parton->SetVeto(ktveto2);
    singlet->push_back(parton);
    parton->SetSing(singlet);
  }
  for (Singlet::const_iterator sit(singlet->begin());
       sit!=singlet->end();++sit) {
    int flow[2]={(*sit)->GetFlow(1),(*sit)->GetFlow(2)};
    if (!(*sit)->GetFlavour().Strong()) continue;
    if ((flow[0]==0 || (*sit)->GetLeft()!=NULL) &&
	(flow[1]==0 || (*sit)->GetRight()!=NULL)) continue;
    if (flow[0]) {
      for (Singlet::const_iterator tit(singlet->begin());
	   tit!=singlet->end();++tit)
	if (tit!=sit && (*tit)->GetFlow(2)==flow[0]) {
	  (*sit)->SetLeft(*tit);
	  (*tit)->SetRight(*sit);
	  break;
	}
    }
    if (flow[1]) {
      for (Singlet::const_iterator tit(singlet->begin());
	   tit!=singlet->end();++tit)
	if (tit!=sit && (*tit)->GetFlow(1)==flow[1]) {
	  (*sit)->SetRight(*tit);
	  (*tit)->SetLeft(*sit);
	  break;
	}
    }
    if ((flow[0] && (*sit)->GetLeft()==NULL) ||
	(flow[1] && (*sit)->GetRight()==NULL))
      THROW(fatal_error,"Missing colour partner");
  }
  for (size_t i(0);i<ampl->Legs().size();++i)
    if (ampl->Leg(i)->K()) {
      Singlet *sing(pmap[ampl->Leg(i)]->GetSing());
      sing->SetSpec(pmap[ampl->IdLeg(ampl->Leg(i)->K())]);
      break;
    }
  return singlet;
}

double CS_Shower::HardScale(const Cluster_Amplitude *const ampl)
{
  if (ampl->Next()) {
    Cluster_Amplitude *next(ampl->Next());
    if (next->NLO()&8) return HardScale(next);
    if (next->OrderQCD()<ampl->OrderQCD()) return ampl->KT2();
    return HardScale(next);
  }
  return ampl->Q2();
}

double CS_Shower::CplFac(const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj,
                         const ATOOLS::Flavour &flk,const int type,
			 const int cpl,const double &mu2) const
{
  cstp::code stp((type&1)?
		 (type&2)?cstp::II:cstp::IF:
		 (type&2)?cstp::FI:cstp::FF);
  return p_shower->GetSudakov()->CplFac(fli, flj, flk,stp,cpl,mu2);
}

double CS_Shower::Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
		       const ATOOLS::Vec4D &pk,const ATOOLS::Flavour &fi,
		       const ATOOLS::Flavour &fj) const
{
  Vec4D npi(pi), npj(pj);
  if (npi[0]<0.0) npi=-pi-pj;
  if (npj[0]<0.0) npj=-pj-pi;
  double pipj(dabs(npi*npj)), pipk(dabs(npi*pk)), pjpk(dabs(npj*pk));
  double Cij(pipk/(pipj+pjpk));
  double Cji(pjpk/(pipj+pipk));
  return 2.0*dabs(pi*pj)/(Cij+Cji);
}

bool CS_Shower::JetVeto(ATOOLS::Cluster_Amplitude *const ampl,
			const int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  msg_Debugging()<<*ampl<<"\n";
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  double q2cut(jf->Ycut()*sqr(rpa->gen.Ecms()));
  size_t noem(0), nospec(0);
  for (size_t i(0);i<ampl->Decays().size();++i) {
    noem|=ampl->Decays()[i]->m_id;
    if (!ampl->Decays()[i]->m_fl.Strong())
      nospec|=ampl->Decays()[i]->m_id;
  }
  msg_Debugging()<<"noem = "<<ID(noem)<<", nospec = "<<ID(nospec)<<"\n";
  double q2min(std::numeric_limits<double>::max());
  size_t imin(0), jmin(0), kmin(0);
  Flavour mofl;
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    if (li->Id()&noem) continue;
    Flavour fi(i<ampl->NIn()?li->Flav().Bar():li->Flav());
    for (size_t j(Max(i+1,ampl->NIn()));j<ampl->Legs().size();++j) {
      Cluster_Leg *lj(ampl->Leg(j));
      if (lj->Id()&noem) continue;
      Flavour fj(j<ampl->NIn()?lj->Flav().Bar():lj->Flav());
      for (size_t k(0);k<ampl->Legs().size();++k) {
	if (k==i || k==j) continue;
	Cluster_Leg *lk(ampl->Leg(k));
	if (lk->Id()&nospec) continue;
	Flavour fk(k<ampl->NIn()?lk->Flav().Bar():lk->Flav());
	cstp::code et((i<ampl->NIn()||j<ampl->NIn())?
		      (k<ampl->NIn()?cstp::II:cstp::IF):
		      (k<ampl->NIn()?cstp::FI:cstp::FF));
	if ((mode==0 && lk->Flav().Strong() &&
	     li->Flav().Strong() && lj->Flav().Strong()) ||
	    p_shower->GetSudakov()->HasKernel(fi,fj,fk,et)) {
	  double q2ijk(Qij2(li->Mom(),lj->Mom(),lk->Mom(),
			    li->Flav(),lj->Flav()));
 	  msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			 <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	  if (q2ijk<0.0) continue;
	  if (mode==0) {
	    if (q2ijk<q2cut) return false;
	  }
	  else {
	    if (q2ijk<q2min) {
	      q2min=q2ijk;
	      mofl=Flavour(kf_gluon);
	      if (li->Flav().IsGluon()) mofl=lj->Flav();
	      if (lj->Flav().IsGluon()) mofl=li->Flav();
	      imin=i;
	      jmin=j;
	      kmin=k;
	    }
	  }
	}
	else {
	  msg_Debugging()<<"No kernel for "<<fi<<" "<<fj
			 <<" <-> "<<fk<<" ("<<et<<")\n";
	}
      }
    }
  }
  if (mode!=0 && imin!=jmin) {
    Vec4D_Vector p=p_cluster->Combine(*ampl,imin,jmin,kmin,mofl,ampl->MS(),1);
    if (p.empty()) {
      msg_Error()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
      return JetVeto(ampl,0);
    }
    Cluster_Amplitude *bampl(Cluster_Amplitude::New());
    bampl->SetProc(ampl->Proc<void>());
    bampl->SetNIn(ampl->NIn());
    bampl->SetJF(ampl->JF<void>());
    for (int i(0), j(0);i<ampl->Legs().size();++i) {
      if (i==jmin) continue;
      if (i==imin) {
	bampl->CreateLeg(p[j],mofl,ampl->Leg(i)->Col());
	bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
      }
      else {
	bampl->CreateLeg(p[j],ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
      }
      ++j;
    }
    bool res=JetVeto(bampl,0);
    bampl->Delete();
    return res;
  }
  msg_Debugging()<<"--- Jet veto ---\n";
  return true;
}

DECLARE_GETTER(CS_Shower,"CSS",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,CS_Shower>::
operator()(const Shower_Key &key) const
{
  return new CS_Shower(key.p_isr,key.p_model,key.p_read);
}

void Getter<Shower_Base,Shower_Key,CS_Shower>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The CSS shower"; 
}

namespace CSSHOWER {

  class CSS_Jet_Criterion: public Jet_Criterion {
  private:
    CS_Shower *p_css;
  public:
    CSS_Jet_Criterion(Shower_Base *const css)
    {
      p_css=dynamic_cast<CS_Shower*>(css);
      if (p_css==NULL) THROW(fatal_error,"CS shower needed but not used");
    }
    bool Jets(Cluster_Amplitude *ampl,int mode)
    {
      return p_css->JetVeto(ampl,mode);
    }
  };// end of class CSS_Jet_Criterion

}// end of namespace CSSHOWER

DECLARE_GETTER(CSS_Jet_Criterion,"CSS",Jet_Criterion,JetCriterion_Key);

Jet_Criterion *Getter<Jet_Criterion,JetCriterion_Key,CSS_Jet_Criterion>::
operator()(const JetCriterion_Key &args) const
{
  return new CSS_Jet_Criterion(args.p_shower);
}

void Getter<Jet_Criterion,JetCriterion_Key,CSS_Jet_Criterion>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The CSS jet criterion";
}
