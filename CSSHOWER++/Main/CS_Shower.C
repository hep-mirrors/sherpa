#include "CSSHOWER++/Main/CS_Shower.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PDF/Main/Jet_Criterion.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

#include <algorithm>
#include <assert.h>

using namespace CSSHOWER;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

CS_Shower::CS_Shower(PDF::ISR_Handler *const _isr,
		     MODEL::Model_Base *const model,
                     const int type) :
  Shower_Base("CSS"), p_isr(_isr), 
  p_shower(NULL), p_cluster(NULL)
{
  Settings& s = Settings::GetMainSettings();
  rpa->gen.AddCitation
    (1,"The Catani-Seymour subtraction based shower is published under \\cite{Schumann:2007mg}.");

  m_maxem = s["CSS_MAXEM"].Get<size_t>();

  SF_Lorentz::SetKappa(s["DIPOLES"]["KAPPA"].Get<double>());

  m_kmode      = s["CSS_KMODE"].Get<int>();
  m_recocheck  = s["CSS_RECO_CHECK"].Get<int>();
  m_respectq2  = s["CSS_RESPECT_Q2"].Get<bool>();

  const int ckfmode { s["CSS_CKFMODE"].Get<bool>() };
  const int pdfcheck{ s["CSS_PDFCHECK"].Get<bool>() };
  const int _qcd    { s["CSS_QCD_MODE"].Get<bool>() };
  const int _qed    { s["CSS_EW_MODE"].Get<bool>() };

  if (_qed==1) {
    s_kftable[kf_photon]->SetResummed();
  }

  p_shower = new Shower(_isr,_qcd,_qed,type);
  p_next = new All_Singlets();
  p_cluster = new CS_Cluster_Definitions(p_shower,m_kmode,pdfcheck,ckfmode);
}

CS_Shower::~CS_Shower() 
{
  CleanUp();
  if (p_shower)      { delete p_shower; p_shower = NULL; }
  if (p_cluster)     { delete p_cluster; p_cluster = NULL; }
  delete p_next;
}

int CS_Shower::PerformShowers(const size_t &maxem,size_t &nem)
{
  if (!p_shower || !m_on) return 1;
  m_weightsmap.Clear();
  Singlet *ls(NULL);
  for (All_Singlets::const_iterator sit(m_allsinglets.begin());
       sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before shower step\n";
    for (Singlet::const_iterator it((*sit)->begin());it!=(*sit)->end();++it)
      if ((*it)->GetPrev()) {
	if (((*it)->GetPrev()->Stat()&2) &&
	    (*it)->GetPrev()==(*it)->GetPrev()->GetSing()->GetSplit()) continue;
	if (m_respectq2) (*it)->SetStart(Min((*it)->KtStart(),(*it)->GetPrev()->KtStart()));
	else (*it)->SetStart((*it)->GetPrev()->KtStart());
      }
    std::map<int,int> colmap;   
    if (ls && (ls->GetSplit()->Stat()&2)) {
      msg_Debugging()<<"Decay. Set color connections.\n";
      Singlet *sing(*sit);
      sing->SetJF(NULL);
      while (sing->GetLeft()) {
	sing=sing->GetLeft()->GetSing();
	sing->SetJF(NULL);
      }
      Parton *d[2]={ls->GetLeft(),ls->GetRight()};
      for (int j=0;j<2;++j) {
	if (d[1-j]->GetFlow(1) || d[1-j]->GetFlow(2)) continue;
	for (int i=1;i<=2;++i) {
	  if (d[j]->GetFlow(i)==0) continue;
	  int nc=Flow::Counter();
	  colmap[nc]=d[j]->GetFlow(i);
	  d[j]->SetFlow(i,nc);
	  d[1-j]->SetFlow(3-i,nc);
	  if (i==1) d[j]->SetLeft(d[1-j]);
	  else d[j]->SetRight(d[1-j]);
	}
      }
      for (Singlet::iterator it((*sit)->begin());it!=(*sit)->end();++it)
	if (*it!=d[0] && *it!=d[1]) (*it)->SetStart(0.0);
    }
    msg_Debugging()<<**sit;
    size_t pem(nem);
    if (!p_shower->EvolveShower(*sit,maxem,nem)) return 0;
    m_weightsmap["PS"] *= p_shower->WeightsMap().at("PS");
    m_weightsmap["PS_QCUT"] *= p_shower->WeightsMap().at("PS_QCUT");
    m_allsinglets=*p_next;
    if (colmap.size()) {
      msg_Debugging()<<"Decay. Reset color connections.\n";
      for (Singlet::iterator
	     it((*sit)->begin());it!=(*sit)->end();++it) {
	for (int i=1;i<=2;++i)
	  if (colmap.find((*it)->GetFlow(i))!=colmap.end()) {
	    if (!(*it)->GetFlavour().Strong()) {
	      (*it)->SetFlow(i,0);
	      continue;
	    }
	    (*it)->SetFlow(i,colmap[(*it)->GetFlow(i)]);
	    (*it)->UpdateColours();
	  }
      }
    }
    msg_Debugging()<<"after shower step with "<<nem-pem
		   <<" of "<<nem<<" emission(s)\n";
    msg_Debugging()<<**sit<<"\n";
    ls=*sit;
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
    if (!(rampl->Flag()&1)) break;
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
  return PrepareStandardShower(ampl);
}

bool CS_Shower::PrepareStandardShower(Cluster_Amplitude *const ampl)
{
  assert(ampl != NULL);
  CleanUp();
  DEBUG_FUNC("");
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
      k->SetStat(cl->Stat());
      almap[apmap[cl]=k]=cl;
      std::map<size_t,Parton*>::iterator cit(kmap.find(cl->Id()));
      if (cit!=kmap.end()) {
	if (k->GetNext()!=NULL) 
	  THROW(fatal_error,"Invalid tree structure. No Next.");
	k->SetNext(cit->second);
	cit->second->SetPrev(k);
	k->SetKScheme(cit->second->KScheme());
	if (k->KScheme()) k->SetMass2(k->Momentum().Abs2());
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
	  else THROW(fatal_error,"Invalid tree structure. LR exist.");
	  kit->second->SetPrev(k);
	  kmap.erase(kit);
	  kit=kmap.begin();
	  split=k;
	}
      kmap[cl->Id()]=k;
    }
    if (sing->GetSpec()) {
      split->SetOldMomentum(split->Momentum());
      sing->SetSpec(sing->GetSpec()->GetNext());
      if (split==NULL) THROW(fatal_error,"Invalid tree structure. No Split.");
      sing->SetSplit(split);
      Parton *l(sing->GetLeft()), *r(sing->GetRight()), *s(sing->GetSpec());
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      split->SetKin(campl->Kin());
      split->SetKScheme((almap[split]->Stat()&4)?1:0);
      if (split->KScheme()) split->SetMass2(split->Momentum().Abs2());
      CS_Parameters cp(p_cluster->KT2
		       (campl->Prev(),almap[l],almap[r],almap[s],
			split->GetType()==pst::FS?split->GetFlavour():
			split->GetFlavour().Bar(),p_ms,
			split->Kin(),split->KScheme(),1));
      cp.m_lt.Invert();
      l->SetLT(cp.m_lt);
      l->SetTest(cp.m_kt2,cp.m_z,cp.m_y,cp.m_phi);
      if (split->KScheme()) split->SetFixSpec(cp.m_pk);
      msg_Debugging()<<"Set reco params: kt = "<<sqrt(cp.m_kt2)<<", z = "
		     <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi
		     <<", mode = "<<cp.m_mode<<", scheme = "<<split->Kin()
		     <<", kmode = "<<split->KScheme()<<"\n";
      sing->SetAll(p_next);
      Vec4D oldl(l->Momentum()), oldr(r->Momentum()), olds(s->Momentum());
      if ((m_recocheck&1) && dynamic_cast<CS_Cluster_Definitions*>
	  (campl->CA<Cluster_Definitions_Base>())!=NULL) {
      std::cout.precision(12);
      Vec4D oldfl(l->FixSpec()), oldfr(r->FixSpec()), oldfs(s->FixSpec());
      Vec4D oldsf(split->FixSpec()), oldso(split->OldMomentum());
      double jcv(0.0);
      p_shower->ReconstructDaughters(sing,jcv,NULL,NULL);
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      CS_Parameters ncp(p_cluster->KT2
			(campl->Prev(),almap[l],almap[r],almap[s],
			 split->GetType()==pst::FS?split->GetFlavour():
			 split->GetFlavour().Bar(),p_ms,
			 split->Kin(),split->KScheme(),1));
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
	if (m_recocheck&2) Abort();
      }
      l->SetMomentum(oldl);
      r->SetMomentum(oldr);
      s->SetMomentum(olds);
      l->SetFixSpec(oldfl);
      r->SetFixSpec(oldfr);
      s->SetFixSpec(oldfs);
      split->SetFixSpec(oldsf);
      split->SetOldMomentum(oldso);
      l->SetLT(cp.m_lt);
      }
      l->SetOldMomentum(oldl);
      r->SetOldMomentum(oldr);
      s->SetOldMomentum(olds);
    }
    double kt2next(0.0);
    if (campl->Prev()) {
      kt2next=campl->Prev()->KT2();
      if (almap[split]->Stat()&2) kt2next=0.0;
    }
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
  p_shower->SetMS(p_ms);
  return true;
}

Singlet *CS_Shower::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap,
 const KT2X_Map &kt2xmap)
{
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  double ktveto2(jf?sqr(jf->Qcut()):
		 sqrt(std::numeric_limits<double>::max()));
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
  singlet->SetJF(jf);
  singlet->SetShower(p_shower);
  singlet->SetMuR2(ampl->MuR2());
  if (ampl->NLO()&8) {
    if (ampl->Next() && 
	(ampl->NIn()+ampl->Leg(2)->NMax()-1>
	 ampl->Legs().size())) ampl->SetNLO(ampl->NLO()&~8);
  }
  singlet->SetNLO(ampl->NLO()&~1);
  if (ampl->Proc<PHASIC::Process_Base>())
    singlet->SetNME(ampl->First()->Legs().size()-ampl->NIn()-
		    ampl->Proc<PHASIC::Process_Base>()->
		    Info().m_fi.NMinExternal());
  if (ampl->Flag()&2) singlet->SetNSkip(1);
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
    parton->SetKin(p_shower->KinScheme());
    if (is) parton->SetBeam(cl->Mom()[3]<0.0?0:1);
    KT2X_Map::const_iterator xit(kt2xmap.find(cl->Id()));
    parton->SetStart(m_respectq2?ampl->MuQ2():xit->second.second);
    if (m_respectq2)
      if (IsDecay(ampl,cl)) parton->SetStart(xit->second.second);
    if (cl->KT2(0)>=0.0) parton->SetSoft(0,cl->KT2(0)); 
    if (cl->KT2(1)>=0.0) parton->SetSoft(1,cl->KT2(1)); 
    if (xit->second.first) singlet->SetNMax(1);
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
	(flow[1] && (*sit)->GetRight()==NULL)) {
      msg_Out()<<METHOD<<" has a problem with\n"<<(*ampl)<<"\n";
      THROW(fatal_error,"Missing colour partner");
    }
  }
  for (size_t i(0);i<ampl->Legs().size();++i)
    if (ampl->Leg(i)->K()) {
      Singlet *sing(pmap[ampl->Leg(i)]->GetSing());
      sing->SetSpec(pmap[ampl->IdLeg(ampl->Leg(i)->K())]);
      break;
    }
  return singlet;
}

int CS_Shower::IsDecay(Cluster_Amplitude *const ampl,
		       Cluster_Leg *const cl) const
{
  for (Cluster_Amplitude *campl(ampl);
       campl;campl=campl->Next()) {
    for (size_t i(0);i<campl->Legs().size();++i) {
      if (campl->Leg(i)->Id()&cl->Id() &&
	  campl->Leg(i)->Stat()&2) return true;
    }
  }
  return false;
}

double CS_Shower::HardScale(const Cluster_Amplitude *const ampl)
{
  if (ampl->Next()) {
    Cluster_Amplitude *next(ampl->Next());
    if (next->NLO()&(4|8)) return HardScale(next);
    if (next->OrderQCD()<ampl->OrderQCD()) return ampl->KT2();
    return HardScale(next);
  }
  return ampl->MuQ2();
}

double CS_Shower::Qij2(const ATOOLS::Vec4D &pi,const ATOOLS::Vec4D &pj,
		       const ATOOLS::Vec4D &pk,const ATOOLS::Flavour &fi,
		       const ATOOLS::Flavour &fj) const
{
  // arXiv:2002.11114 [hep-ph]
  const double beta(0.5);
  double t1(2.0*(pi*pj)*(pj*pk)/(pi*pk));
  double t2(2.0*(pj*pi)*(pi*pk)/(pj*pk));
  double xi1(dabs((pi*pj)/(pk*pj)));
  double xi2(dabs((pj*pi)/(pk*pi)));
  t1*=pow(Max(xi1,1.0/xi1),-beta/2.0);
  t2*=pow(Max(xi2,1.0/xi2),-beta/2.0);
  if (pi[0]<0.0) return dabs(t1);
  if (pj[0]<0.0) return dabs(t2);
  return Min(t1,t2);
}

double CS_Shower::JetVeto(ATOOLS::Cluster_Amplitude *const ampl,
			  const int mode)
{
  DEBUG_FUNC("mode = "<<mode);
  msg_Debugging()<<*ampl<<"\n";
  NLO_subevtlist *subs(NULL);
  if (mode) subs=ampl->Proc<PHASIC::Process_Base>()->GetRSSubevtList();
  size_t noem(0), nospec(0);
  for (size_t i(0);i<ampl->Decays().size();++i) {
    noem|=ampl->Decays()[i]->m_id;
    if (!ampl->Decays()[i]->m_fl.Strong())
      nospec|=ampl->Decays()[i]->m_id;
  }
  msg_Debugging()<<"noem = "<<ID(noem)<<", nospec = "<<ID(nospec)<<"\n";
  std::set<Q2_Value> q2list;
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
	if (subs) {
	  bool found(false);
	  for (size_t l(0);l<subs->size()-1;++l) {
	    NLO_subevt *sub((*subs)[l]);
	    if (k==sub->m_k && ((i==sub->m_i && j==sub->m_j) ||
				(i==sub->m_j && j==sub->m_i))) {
	      found=true;
	      break;
	    }
	  }
	  if (!found) continue;
	}
	Cluster_Leg *lk(ampl->Leg(k));
	if (lk->Id()&nospec) continue;
	Flavour fk(k<ampl->NIn()?lk->Flav().Bar():lk->Flav());
	if (lk->Flav().Strong() &&
	    li->Flav().Strong() && lj->Flav().Strong() &&
	    (li->Flav().IsGluon() || lj->Flav().IsGluon() ||
	     li->Flav()==lj->Flav().Bar())) {
	  double q2ijk(Qij2(li->Mom(),lj->Mom(),lk->Mom(),
			    li->Flav(),lj->Flav()));
	  msg_Debugging()<<"Q_{"<<ID(li->Id())<<ID(lj->Id())
			 <<","<<ID(lk->Id())<<"} = "<<sqrt(q2ijk)<<"\n";
	  if (q2ijk<0.0) continue;
	  Flavour mofl=Flavour(kf_gluon);
	  if (li->Flav().IsGluon()) mofl=lj->Flav();
	  if (lj->Flav().IsGluon()) mofl=li->Flav();
	  q2list.insert(Q2_Value(q2ijk,mofl,i,j,k));
	}
	else {
	  msg_IODebugging()<<"No kernel for "<<fi<<" "<<fj<<" <-> "<<fk<<"\n";
	}
      }
    }
  }
  if (mode==0) {
    double q2min(std::numeric_limits<double>::max());
    if (q2list.size()) q2min=q2list.begin()->m_q2;
    msg_Debugging()<<"--- "<<sqrt(q2min)<<" ---\n";
    return q2min;
  }
  while (q2list.size()) { 
    Flavour mofl(q2list.begin()->m_fl);
    size_t imin(q2list.begin()->m_i);
    size_t jmin(q2list.begin()->m_j);
    size_t kmin(q2list.begin()->m_k);
    q2list.erase(q2list.begin());
    Cluster_Param cp=p_cluster->Cluster
      (Cluster_Config(ampl,imin,jmin,kmin,mofl,ampl->MS(),NULL,1));
    if (cp.m_pijt==Vec4D())
      cp=p_cluster->Cluster(Cluster_Config(ampl,jmin,imin,kmin,mofl,ampl->MS(),NULL,1));
    if (cp.m_pijt==Vec4D()) continue;
    Cluster_Amplitude *bampl(Cluster_Amplitude::New());
    bampl->SetProc(ampl->Proc<void>());
    bampl->SetNIn(ampl->NIn());
    bampl->SetJF(ampl->JF<void>());
    for (int i(0), j(0);i<ampl->Legs().size();++i) {
      if (i==jmin) continue;
      if (i==imin) {
	bampl->CreateLeg(cp.m_pijt,mofl,ampl->Leg(i)->Col());
	bampl->Legs().back()->SetId(ampl->Leg(imin)->Id()|ampl->Leg(jmin)->Id());
	bampl->Legs().back()->SetK(ampl->Leg(kmin)->Id());	
      }
      else {
	bampl->CreateLeg(i==kmin?cp.m_pkt:cp.m_lam*ampl->Leg(i)->Mom(),
			 ampl->Leg(i)->Flav(),ampl->Leg(i)->Col());
      }
      ++j;
    }
    double res=JetVeto(bampl,0);
    bampl->Delete();
    if (res==std::numeric_limits<double>::max()) continue;
    return res;
  }
  msg_Debugging()<<METHOD<<"(): Combine failed. Use R configuration."<<std::endl;
  return JetVeto(ampl,0);
}

DECLARE_GETTER(CS_Shower,"CSS",Shower_Base,Shower_Key);

Shower_Base *Getter<Shower_Base,Shower_Key,CS_Shower>::
operator()(const Shower_Key &key) const
{
  return new CS_Shower(key.p_isr,key.p_model,key.m_type);
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
    double Value(Cluster_Amplitude *ampl,int mode)
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
