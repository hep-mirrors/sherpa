#ifndef CSSHOWER_Main_CS_Shower_H
#define CSSHOWER_Main_CS_Shower_H

#include "PDF/Main/Shower_Base.H"

#include "PDF/Main/ISR_Handler.H"
#include "CSSHOWER++/Showers/Shower.H"
#include "CSSHOWER++/Tools/Singlet.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Blob.H"

#include "CSSHOWER++/Main/CS_Cluster_Definitions.H"
#include "PHASIC++/Selectors/Jet_Finder.H"

namespace ATOOLS { 
  class Cluster_Leg;
  class Data_Reader; 
}

namespace CSSHOWER {
  class Splitting_Function_Base;
}

namespace CSSHOWER {

  typedef std::map<size_t,std::pair<double,double> > KT2X_Map;

  class CS_Shower : public PDF::Shower_Base {
  
  private:
 
    PDF::ISR_Handler * p_isr;
    int m_weightmode, m_kmode, m_dmode;
    size_t m_maxem;
    
    Shower          * p_shower;
    All_Singlets m_allsinglets;
    CS_Cluster_Definitions *p_cluster;
    All_Singlets *p_next;

    ATOOLS::Mass_Selector *p_ms;

    ATOOLS::Cluster_Amplitude *p_ampl, *p_rampl;
    ATOOLS::Particle_List      m_psp;

    size_t IdCount(const size_t &id);
    void   GetKT2Min(ATOOLS::Cluster_Amplitude *const ampl,const size_t &id,
		     KT2X_Map &kt2xmap,std::set<size_t> &aset);
    void   GetKT2Min(ATOOLS::Cluster_Amplitude *const ampl,KT2X_Map &kt2xmap);

    double HardScale(const ATOOLS::Cluster_Amplitude *const ampl);

    double CouplingWeight(const size_t &oqcd,const double &kt2,
			  const double &kt2r) const;
    
    Singlet *TranslateAmplitude(ATOOLS::Cluster_Amplitude *const ampl,
				std::map<ATOOLS::Cluster_Leg*,Parton*> &pmap,
				std::map<Parton*,ATOOLS::Cluster_Leg*> &lmap,
				const KT2X_Map &kt2xmap);

    int PerformShowers(const size_t &maxem,size_t &nem);

  public:

    // constructor 
    CS_Shower(PDF::ISR_Handler *const isr, MODEL::Model_Base *const model,
	      ATOOLS::Data_Reader *const dataread); 

    // destructor
    ~CS_Shower();

    //member functions
    int  PerformShowers();
    int  PerformDecayShowers();
    int  TrialEmission();

    double CouplingWeight(ATOOLS::Cluster_Amplitude *const ampl);
    double TrialWeight(ATOOLS::Cluster_Amplitude *const ampl);

    bool ExtractPartons(ATOOLS::Blob_List *const blist);

    void CleanUp();

    // inline functions
    ATOOLS::Cluster_Definitions_Base * GetClusterDefinitions();
    ATOOLS::Cluster_Amplitude *GetRealEmissionAmplitude();
    bool PrepareShower(ATOOLS::Cluster_Amplitude *const ampl);
    double CalculateWeight(ATOOLS::Cluster_Amplitude *const ampl);      
    double CalculateAnalyticWeight(ATOOLS::Cluster_Amplitude *const ampl);      

    std::string GetKT2(const std::string &jm2) const;

  };
}

#endif

#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

CS_Shower::CS_Shower(PDF::ISR_Handler *const _isr,MODEL::Model_Base *const model,
		     Data_Reader *const _dataread) : 
  Shower_Base("CSS"), p_isr(_isr), 
  p_shower(NULL), p_cluster(NULL), p_ampl(NULL)
{
  rpa.gen.AddCitation
    (1,"The Catani-Seymour subtraction based shower is published under \\cite{Schumann:2007mg}.");
  int maxem=_dataread->GetValue<int>("CSS_MAXEM",-1);
  if (maxem<0) m_maxem=std::numeric_limits<size_t>::max();
  else {
    m_maxem=maxem;
    msg_Info()<<METHOD<<"(): Set max emissions "<<m_maxem<<"\n";
  }
  m_kmode=_dataread->GetValue<int>("CSS_KMODE",1);
  if (m_kmode!=1) msg_Info()<<METHOD<<"(): Set kernel mode "<<m_kmode<<"\n";
  m_dmode=_dataread->GetValue<int>("CSS_DMODE",0);
  if (m_dmode!=0) msg_Info()<<METHOD<<"(): Set DIS mode "<<m_dmode<<"\n";
  
  m_weightmode = int(_dataread->GetValue<int>("WEIGHT_MODE",1));
  
  int _qed=_dataread->GetValue<int>("CSS_EW_MODE",3);
  p_shower = new Shower(_isr,_qed,_dataread);
  
  p_next = new All_Singlets();
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
  if (!p_shower) return 1;
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    msg_Debugging()<<"before shower step\n";
      for (Singlet::const_iterator it((*sit)->begin());it!=(*sit)->end();++it)
	if ((*it)->GetPrev()) 
	  if((*it)->GetPrev()->GetNext()==*it)
	    (*it)->SetStart((*it)->GetPrev()->KtStart());
      msg_Debugging()<<**sit;
    if (!p_shower->EvolveShower(*sit,maxem,nem)) return 0;
    if ((*sit)->GetLeft()) p_shower->ReconstructDaughters(*sit,true);
    msg_Debugging()<<"after shower step\n";
      msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  return 1;
}

int CS_Shower::PerformShowers() 
{
  size_t nem(0);
  return PerformShowers(m_maxem,nem);
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
  
  for (size_t i(0);i<m_psp.size();++i)
    psblob->AddToInParticles(m_psp[i]);
  m_psp.clear();

  psblob->SetStatus(blob_status::needs_beams |
		    blob_status::needs_harddecays |
		    blob_status::needs_hadronization);
  
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit)
      (*sit)->ExtractPartons(psblob,p_ms);
  return true;
}

void CS_Shower::CleanUp()
{
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
    if (*sit) delete *sit;
  }
  m_allsinglets.clear();
  for (size_t i(0);i<m_psp.size();++i)
    delete m_psp[i];
  m_psp.clear();
}

ATOOLS::Cluster_Definitions_Base * CS_Shower::GetClusterDefinitions() 
{
  if (p_cluster==NULL)
    p_cluster = new CS_Cluster_Definitions(p_shower,m_kmode);
  return p_cluster;
}

size_t CS_Shower::IdCount(const size_t &id)
{
  size_t ic(id), cn(0);
  for (size_t i(0);ic>0;++i) {
    size_t c(1<<i);
    if (ic&c) {
      ++cn;
      ic-=c;
    }
  }
  return cn;
}

void CS_Shower::GetKT2Min(Cluster_Amplitude *const ampl,const size_t &id,
			  KT2X_Map &kt2xmap,std::set<size_t> &aset)
{
  msg_Indent();
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cl(ampl->Leg(i));
    if ((cl->Id()&id)==0) continue;
    if (IdCount(cl->Id())>1) GetKT2Min(ampl->Prev(),cl->Id(),kt2xmap,aset);
    else kt2xmap[cl->Id()].first=kt2xmap[cl->Id()].second=HardScale(ampl);
    if (cl->K()) {
      std::vector<size_t> cns;
      double ckt2min(std::numeric_limits<double>::max()), ckt2max(0.0);
      kt2xmap[cl->Id()].first=kt2xmap[cl->Id()].second=HardScale(ampl->Prev());
      for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
	if (kit->first!=cl->Id() && kit->first&cl->Id() &&
	    aset.find(kit->first)==aset.end()) {
	  ckt2min=Min(ckt2min,kit->second.first);
	  ckt2max=Max(ckt2max,kit->second.second);
	  if (cl->Stat()==3) {
	    bool ins(true);
	    for (size_t j(0);j<cns.size();++j)
	      if (cns[j]&kit->first) {
		ins=false;
		break;
	      }
	    if (ins) cns.push_back(kit->first);
	  }
	}
      bool smin(true);
      if (cl->Stat()==3 && cns.size()<cl->DMax()) smin=false;
      for (KT2X_Map::iterator kit(kt2xmap.begin());kit!=kt2xmap.end();++kit)
	if (kit->first!=cl->Id() && kit->first&cl->Id() &&
	    aset.find(kit->first)==aset.end()) {
	  if (smin) kit->second.first=ckt2min;
	  else kit->second.first=0.0;
	  kit->second.second=ckt2max;
	  if (cl->Stat()==3) aset.insert(kit->first);
	}
    }
    if (cl->Stat()==3) {
      kt2xmap[cl->Id()].first=std::numeric_limits<double>::max();
      kt2xmap[cl->Id()].second=0.0;
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
  bool smin(cns.size()-ampl->NIn()==campl->Leg(2)->NMax());
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

bool CS_Shower::PrepareShower(Cluster_Amplitude *const ampl)
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
      sing->SetSpec(sing->GetSpec()->GetNext());
      if (split==NULL) THROW(fatal_error,"Invalid tree structure");
      sing->SetSplit(split);
      Parton *l(sing->GetLeft()), *r(sing->GetRight()), *s(sing->GetSpec());
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      CS_Parameters cp(p_cluster->KT2
		       (almap[l],almap[r],almap[s],
			split->GetFlavour(),p_ms));
      l->SetTest(cp.m_kt2,cp.m_z,cp.m_y,cp.m_phi);
      l->SetStart(cp.m_kt2);
      r->SetStart(cp.m_kt2);
      msg_Debugging()<<"Set reco params: kt = "<<sqrt(cp.m_kt2)<<", z = "
		     <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi
		     <<", mode = "<<cp.m_mode<<"\n";
      sing->SetAll(p_next);
#ifdef USING__reco_check
      std::cout.precision(12);
      Vec4D oldl(l->Momentum()), oldr(r->Momentum()), olds(s->Momentum());
      sing->BoostBackAllFS(l,r,s,split,split->GetFlavour(),cp.m_mode);
      p_shower->ReconstructDaughters(sing,true);
      almap[l]->SetMom(almap[l]->Id()&3?-l->Momentum():l->Momentum());
      almap[r]->SetMom(almap[r]->Id()&3?-r->Momentum():r->Momentum());
      almap[s]->SetMom(almap[s]->Id()&3?-s->Momentum():s->Momentum());
      CS_Parameters ncp(p_cluster->KT2
			(almap[l],almap[r],almap[s],
			 split->GetFlavour(),p_ms));
      msg_Debugging()<<"New reco params: kt = "<<sqrt(ncp.m_kt2)<<", z = "
		     <<ncp.m_z<<", y = "<<ncp.m_y<<", phi = "<<ncp.m_phi<<"\n";
      msg_Debugging()<<"            vs.: kt = "<<sqrt(cp.m_kt2)<<", z = "
		     <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi<<"\n";
      if (!IsEqual(ncp.m_kt2,cp.m_kt2,1.0e-6) || 
	  !IsEqual(ncp.m_z,cp.m_z,1.0e-6) || 
	  !IsEqual(ncp.m_y,cp.m_y,1.0e-6) || 
	  !IsEqual(ncp.m_phi,cp.m_phi,1.0e-6) ||
	  !IsEqual(oldl,l->Momentum(),1.0e-6) || 
	  !IsEqual(oldr,r->Momentum(),1.0e-6) || 
	  !IsEqual(olds,s->Momentum(),1.0e-6)) {
	msg_Error()<<"\nFaulty reco params: kt = "<<sqrt(ncp.m_kt2)<<", z = "
		   <<ncp.m_z<<", y = "<<ncp.m_y<<", phi = "<<ncp.m_phi<<"\n";
	msg_Error()<<"               vs.: kt = "<<sqrt(cp.m_kt2)<<", z = "
		   <<cp.m_z<<", y = "<<cp.m_y<<", phi = "<<cp.m_phi<<"\n";
	msg_Error()<<"  "<<oldl<<" "<<oldr<<" "<<olds<<"\n";
	msg_Error()<<"  "<<l->Momentum()<<" "<<r->Momentum()
		   <<" "<<s->Momentum()<<"\n";
      }
      l->SetMomentum(oldl);
      r->SetMomentum(oldr);
      s->SetMomentum(olds);
#endif
      sing->BoostBackAllFS(l,r,s,split,split->GetFlavour(),cp.m_mode);
    }
    double kt2prev(campl->Next()?campl->KT2QCD():kt2xmap[1].second);
    double kt2next(campl->Prev()?campl->Prev()->KT2QCD():0.0);
    for (size_t i(0);i<campl->Legs().size();++i) {
      std::map<Cluster_Leg*,Parton*>::const_iterator 
	pit(apmap.find(campl->Leg(i)));
      if (pit!=apmap.end()) {
	pit->second->SetKtPrev(kt2prev);
	pit->second->SetKtNext(kt2next);
      }
    }
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
  for (All_Singlets::const_iterator 
	 sit(m_allsinglets.begin());sit!=m_allsinglets.end();++sit) {
      for (Singlet::const_iterator 
	     pit((*sit)->begin());pit!=(*sit)->end();++pit) {
	if ((*pit)->GetPrev()) {
	  if ((*pit)->GetPrev()->GetNext()==*pit) 
	    (*pit)->SetStart((*pit)->GetPrev()->KtStart());
	  (*pit)->SetKtPrev((*pit)->GetPrev()->KtNext());
	}
      }
      (*sit)->SetJF(ampl->JF<PHASIC::Jet_Finder>());
      (*sit)->SetAll(p_next);
      msg_Debugging()<<**sit;
    msg_Debugging()<<"\n";
  }
  msg_Debugging()<<"}\n";
  p_shower->SetRBMax(ampl->RBMax());
  p_shower->SetMS(p_ms);
  return true;
}

Singlet *CS_Shower::TranslateAmplitude
(Cluster_Amplitude *const ampl,
 std::map<Cluster_Leg*,Parton*> &pmap,std::map<Parton*,Cluster_Leg*> &lmap,
 const KT2X_Map &kt2xmap)
{
  PHASIC::Jet_Finder *jf(ampl->JF<PHASIC::Jet_Finder>());
  double ktveto2(jf?jf->GlobalYcut()*sqr(rpa.gen.Ecms()):4.0*ampl->MuR2());
  Singlet *singlet(new Singlet());
  singlet->SetMS(p_ms);
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
    pmap[cl]=parton;
    lmap[parton]=cl;
    parton->SetRFlow();
    if (is) {
      if (Vec3D(p.Momentum())*Vec3D(rpa.gen.PBeam(0))>0.) {
	parton->SetXbj(p.Momentum()[0]/rpa.gen.PBeam(0)[0]);
	parton->SetBeam(0);
      }
      else { 
	parton->SetXbj(p.Momentum()[0]/rpa.gen.PBeam(1)[0]);
	parton->SetBeam(1);
      }
    }
    KT2X_Map::const_iterator xit(kt2xmap.find(cl->Id()));
    if (cl->Q2Shower()<0.0) parton->SetStart(xit->second.second);
    else parton->SetStart(cl->Q2Shower());
    parton->SetKtMax(xit->second.first);
    parton->SetVeto(ktveto2);
    singlet->push_back(parton);
    parton->SetSing(singlet);
  }
  for (Singlet::const_iterator sit(singlet->begin());
       sit!=singlet->end();++sit) {
    if (!(*sit)->GetFlavour().Strong() ||
	(*sit)->GetLeft() || (*sit)->GetRight()) continue;
    int flow[2]={(*sit)->GetFlow(1),(*sit)->GetFlow(2)};
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
	(flow[0] && (*sit)->GetLeft()==NULL))
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


double CS_Shower::CalculateWeight(Cluster_Amplitude *const ampl)
{
  if (ampl->JF<PHASIC::Jet_Finder>()==NULL) return 1.0;
  switch (m_weightmode) {
  case 1 : return CalculateAnalyticWeight(ampl);
  default : THROW(fatal_error,"Invalid reweighting mode"); 
  }
  return 1.0;
}


double CS_Shower::CalculateAnalyticWeight(Cluster_Amplitude *const ampl) 
{
  // calculate weight
  p_ms=ampl->MS();
  msg_Debugging()<<METHOD<<"(): {\n";
  double wgt(1.0);
  {
    msg_Indent();
    std::map<size_t,Cluster_Leg*> legs;
    // find hard scale
    double kt2max(ampl->MuF2());
    Cluster_Amplitude *ref(ampl);
    while (ref->Next()) {
      ref=ref->Next();
      for (size_t i(0);i<ref->Legs().size();++i) {
	const Cluster_Leg *cl(ref->Leg(i));  
	//to be changed!!!!
	kt2max=Max(kt2max,dabs(cl->Mom().Abs2()));
      }
    }
    kt2max=Max(kt2max,HardScale(ref));
    msg_Debugging()<<"largest ps scale "<<sqrt(kt2max)<<"\n";
    for (size_t i(0);i<ref->Legs().size();++i) {
      legs[ref->Leg(i)->Id()]=ref->Leg(i);
    }
    if (ref->OrderQCD()) {
      double cf((ref->Leg(0)->Flav().Strong()||
		 ref->Leg(ref->NIn()-1)->Flav().Strong())?
		p_shower->GetSudakov()->ISCplFac():
		p_shower->GetSudakov()->FSCplFac());
      msg_Debugging()<<"core => mu = "<<sqrt(cf)
		     <<"*"<<sqrt(ref->MuF2())<<" {\n";
      {
	msg_Indent();
	wgt*=CouplingWeight(ref->OrderQCD(),cf*ref->MuF2(),ref->MuR2());
      }
      msg_Debugging()<<"}\n";
    }
    size_t istag(((1<<ampl->NIn())-1));
    while (ref->Prev()) {
      ref=ref->Prev();
      bool split(false);
      for (size_t i(0);i<ref->Legs().size();++i) {
	size_t idi(ref->Leg(i)->Id());
	if (legs.find(idi)==legs.end()) {
	  for (size_t j(i+1);j<ref->Legs().size();++j) {
	    size_t idj(ref->Leg(j)->Id());
	    if (legs.find(idi+idj)!=legs.end()) {
	      double cf((idi&istag)?
			p_shower->GetSudakov()->ISCplFac():
			p_shower->GetSudakov()->FSCplFac());
	      double ckt2(ref->KT2QCD());
	      const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFMap());
	      if ((idi+idj)&istag) {
		if (legs[idi+idj]->K()&istag) cmap=&p_shower->GetSudakov()->IIMap();
		else cmap=&p_shower->GetSudakov()->IFMap();
	      }
	      else {
		if (legs[idi+idj]->K()&istag) cmap=&p_shower->GetSudakov()->FIMap();
	      }
	      SF_EEE_Map::const_iterator eees(cmap->find(ref->Leg(i)->Flav()));
	      if (eees!=cmap->end()) {
		SF_EE_Map::const_iterator ees(eees->second.find(ref->Leg(j)->Flav()));
		if (ees!=eees->second.end()) {
		  SF_E_Map::const_iterator es(ees->second.find(legs[idi+idj]->Flav()));
		  if (es!=ees->second.end()) {
		    cf=es->second->Coupling()->CplFac(ckt2);
		    msg_Debugging()<<"known ";
		  }
		}
	      }
	      msg_Debugging()<<"split "<<ID(idi+idj)
			<<" -> "<<ID(idi)<<","<<ID(idj)
			<<" => kt = "<<sqrt(cf)<<"*"<<sqrt(ckt2)<<" {\n";
	      msg_Indent();
	      size_t coqcd(ref->OrderQCD()-ref->Next()->OrderQCD());
	      if (coqcd) wgt*=CouplingWeight(coqcd,cf*ckt2,ref->MuR2());
	      legs[idi]=ref->Leg(i);
	      legs[idj]=ref->Leg(j);
	      split=true;
	      break;
	      
	    }
	  }
	  msg_Debugging()<<"}\n";
	  break;
	}
      }
      if (!split) THROW(fatal_error,"Internal error");
    } 
  }
  msg_Debugging()<<"} -> w = "<<wgt<<"\n";
  return wgt;
}

double CS_Shower::HardScale(const Cluster_Amplitude *const ampl)
{
  if (ampl->Next()) {
    Cluster_Amplitude *next(ampl->Next());
    if (next->OrderQCD()<ampl->OrderQCD()) return ampl->KT2QCD();
    for (size_t i(0);i<next->Legs().size();++i)
      if (next->Leg(i)->K()) {
	Vec4D sum;
	size_t cid(next->Leg(i)->Id());
	for (size_t j(0);j<ampl->Legs().size();++j)
	  if (ampl->Leg(j)->Id()&cid) sum+=ampl->Leg(j)->Mom();
	return dabs(sum.Abs2());
      }
    // this should catch the partonic hadron decays
    if (ampl->NIn()==1) return (ampl->Leg(0)->Mom()).Abs2();
    THROW(fatal_error,"Invalid amplitude");
  }
  double q2cut(0.0);
  if (ampl->JF<PHASIC::Jet_Finder>()) {
    q2cut=ampl->JF<PHASIC::Jet_Finder>()->GlobalYcut();
    q2cut*=sqr(rpa.gen.Ecms());
  }
  double xf(m_dmode?ampl->X1()*ampl->X2():1.0);
  return Max(ampl->MuF2()/xf,q2cut);
}

double CS_Shower::CouplingWeight(const size_t &oqcd,const double &kt2,
				 const double &kt2r) const
{
  double asc((*MODEL::as)(kt2));
  double asr((*MODEL::as)(kt2r));
  msg_Debugging()<<"as weight (\\alpha_s("<<sqrt(kt2)
		 <<")/\\alpha_s("<<sqrt(kt2r)<<"))^O_{as} = ( "
		 <<asc<<" / "<<asr<<" ) ^ "<<oqcd<<" = "
		 <<pow(asc/asr,int(oqcd))<<"\n";
  return pow(asc/asr,int(oqcd));
}

ATOOLS::Cluster_Amplitude *CS_Shower::GetRealEmissionAmplitude()
{
  if (p_ampl) p_ampl->Delete();
  Cluster_Amplitude *ampl(p_rampl);
  while (ampl->Next()) ampl=ampl->Next();
  p_ampl=ampl->Copy();
  ampl=p_ampl->InitNext();
  ampl->CopyFrom(p_ampl);
  Parton *const *split(p_shower->GetLast());
  for (int i(0);i<3;++i)
    if (split[i]->GetType()==pst::FS) {
      if (split[i]->Id()==0) {
	ampl->CreateLeg(split[i]->Momentum(),split[i]->GetFlavour(),
			ColorID(split[i]->GetFlow(1),split[i]->GetFlow(2)));
	ampl->Legs().back()->SetQ2Shower(split[i]->KtStart());
	ampl->Legs().back()->SetStat(1);
      }
      else {
	Cluster_Leg *cl(ampl->IdLeg(split[i]->Id()));
	cl->SetMom(split[i]->Momentum());
	cl->SetFlav(split[i]->GetFlavour());
	cl->SetCol(ColorID(split[i]->GetFlow(1),split[i]->GetFlow(2)));
	cl->SetQ2Shower(split[i]->KtStart());
	p_ampl->IdLeg(cl->Id())->SetQ2Shower(split[i]->KtStart());
      }
    }
    else {
      if (split[i]->Id()==0) {
	ampl->CreateLeg(-split[i]->Momentum(),split[i]->GetFlavour().Bar(),
			ColorID(split[i]->GetFlow(2),split[i]->GetFlow(1)));
	ampl->Legs().back()->SetQ2Shower(split[i]->KtStart());
	ampl->Legs().back()->SetStat(1);
      }
      else {
	Cluster_Leg *cl(ampl->IdLeg(split[i]->Id()));
	cl->SetMom(-split[i]->Momentum());
	cl->SetFlav(-split[i]->GetFlavour());
	cl->SetCol(ColorID(split[i]->GetFlow(2),split[i]->GetFlow(1)));
	cl->SetQ2Shower(split[i]->KtStart());
	p_ampl->IdLeg(cl->Id())->SetQ2Shower(split[i]->KtStart());
      }
    }
  ampl->SetIdNew(1<<p_ampl->Legs().size());
  return p_ampl;
}

int CS_Shower::TrialEmission()
{
  size_t nem(0);
  int res(PerformShowers(1,nem));
  if (res!=1) return res;
  return nem==1;
}

double CS_Shower::TrialWeight(ATOOLS::Cluster_Amplitude *const ampl)
{
  double as((*MODEL::as)(rpa.gen.CplScale()));
  if (p_cluster==NULL) p_cluster = new CS_Cluster_Definitions(p_shower,m_kmode);
  DEBUG_FUNC("");
  p_ms=ampl->MS();
  p_shower->SetMS(p_ms);
  Cluster_Amplitude *next(ampl->Next());
  Cluster_Leg *lj(next->IdLeg(next->IdNew()));
  msg_Debugging()<<*ampl<<"\n"<<*next<<"\nnew = "<<*lj<<"\n";
  double wgt(0.0);
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *lijt(ampl->Leg(i)), *li(next->IdLeg(lijt->Id())); 
    if (!lijt->Flav().Strong()) continue;
    for (size_t k(0);k<ampl->Legs().size();++k) {
      if (i==k) continue;
      Cluster_Leg *lkt(ampl->Leg(k));
      if (!lkt->Flav().Strong()) continue;
      const SF_EEE_Map *cmap(&p_shower->GetSudakov()->FFMap());
      if (i<ampl->NIn()) {
	if (k<ampl->NIn()) cmap=&p_shower->GetSudakov()->IIMap();
	else cmap=&p_shower->GetSudakov()->IFMap();
      }
      else {
	if (k<ampl->NIn()) cmap=&p_shower->GetSudakov()->FIMap();
	else cmap=&p_shower->GetSudakov()->FFMap();
      }
      SF_EEE_Map::const_iterator eees(cmap->find(li->Flav()));
      if (eees==cmap->end()) {
	msg_Error()<<METHOD<<"(): No SF for emitted flavour "<<li->Flav()<<".\n";
	continue;
      }
      SF_EE_Map::const_iterator ees(eees->second.find(lj->Flav()));
      if (ees==eees->second.end()) {
	msg_Error()<<METHOD<<"(): No SF for emitted flavour "<<lj->Flav()<<".\n";
	continue;
      }
      SF_E_Map::const_iterator es(ees->second.find(lijt->Flav()));
      if (es==ees->second.end()) {
	msg_Debugging()<<METHOD<<"(): No SF for emitter flavour "
		       <<lijt->Flav()<<".\n";
	continue;
      }
      Splitting_Function_Base *cdip(es->second);
      cdip->SetFlavourSpec(k<ampl->NIn()?lkt->Flav().Bar():lkt->Flav());
      double Q2=(lijt->Mom()+lkt->Mom()).Abs2(), scale=Q2;
      Cluster_Leg *lk(next->IdLeg(lkt->Id()));
      CS_Parameters cs(p_cluster->KT2(li,lj,lk,lijt->Flav(),p_ms));
      double eta=1.0;// temporary, FF only
      double w=8.0*M_PI*as/((li->Mom()+lj->Mom()).Abs2()
			    -sqr(p_ms->Mass(lijt->Flav())))
	*(*cdip)(cs.m_z,cs.m_y,eta,scale,Q2,1);
      msg_Debugging()<<"add  {\n  "<<*lijt<<"\n  "<<*lkt
		     <<"\n} -> w = "<<w<<"\n";
      wgt+=w;
    }
  }
  msg_Debugging()<<"weight = "<<wgt<<"\n";
  return wgt*ampl->RBMax();
}

double CS_Shower::CouplingWeight(ATOOLS::Cluster_Amplitude *const ampl)
{
  double kt2(ampl->KT2QCD());
  size_t idi(ampl->Next()->IdLeg(ampl->Next()->IdNew())->Id());
  double cf((idi&((1<<ampl->NIn())-1))?
	    p_shower->GetSudakov()->ISCplFac():
	    p_shower->GetSudakov()->FSCplFac());
  if (kt2/cf<p_shower->GetSudakov()->PT2Min()) return 1.0;
  size_t coqcd(ampl->OrderQCD()-ampl->Next()->OrderQCD());
  if (coqcd) return CouplingWeight(coqcd,cf*kt2,ampl->MuR2());
  return 1.0;
}

std::string CS_Shower::GetKT2(const std::string &jm2) const
{
  return "0.25*"+jm2;
}

namespace PDF {

  DECLARE_GETTER(CSS_Getter,"CSS",Shower_Base,Shower_Key);

  Shower_Base *CSS_Getter::operator()(const Shower_Key &key) const
  {
    return new CS_Shower(key.p_isr,key.p_model,key.p_read);
  }

  void CSS_Getter::PrintInfo(std::ostream &str,const size_t width) const
  { 
    str<<"The CSS shower"; 
  }

}
