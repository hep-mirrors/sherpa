#include "POWHEG/Main/CS_Gamma.H"

#include "POWHEG/Main/CS_POWHEG.H"
#include "POWHEG/Showers/Splitting_Function_Base.H"
#include "PHASIC++/Process/POWHEG_Process.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"

using namespace POWHEG;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

#define DEBUG__Trial_Weight

Weight_Key::Weight_Key(const size_t &ij,const size_t &k,
		       const ATOOLS::Flavour &flij,
		       const ATOOLS::Flavour &fli,
		       const ATOOLS::Flavour &flj):
  m_ij(ij), m_k(k), m_rbkey(((m_ij&1||m_ij&2)|2*(m_k&1||m_k&2)),flij,fli,flj)
{
}

CS_Gamma::CS_Gamma(CS_POWHEG *const css,Shower *const shower,
		   CS_Cluster_Definitions *const cluster):
  p_css(css), p_shower(shower), p_cluster(cluster),
  p_rproc(NULL), m_on(0), m_zhth(1.0e6), m_ktres(0.0)
{
}

Weight_Map CS_Gamma::CalculateWeight(Cluster_Amplitude *const ampl)
{
  Cluster_Amplitude *rampl(ampl->Copy());
#ifdef DEBUG__Trial_Weight
  DEBUG_FUNC(ampl);
  msg_Debugging()<<*rampl<<"\n";
#endif
  rampl->SetIdNew(0);
  std::map<size_t,size_t> idmap;
  for (size_t i(0);i<rampl->Legs().size();++i) {
    idmap[1<<i]=rampl->Leg(i)->Id();
    rampl->Leg(i)->SetId(1<<i);
  }
  for (size_t j(0);j<rampl->Decays().size();++j) {
    size_t oid(rampl->Decays()[j]->m_id), nid(oid);
    for (std::map<size_t,size_t>::const_iterator
	   iit(idmap.begin());iit!=idmap.end();++iit)
      if (oid&iit->second) {
	nid&=~iit->second;
	nid|=iit->first;
      }
    rampl->Decays()[j]->m_id=nid;
  }
  std::map<nlo_type::code,Process_Map*> *procs
    (ampl->Procs<std::map<nlo_type::code,Process_Map*> >());
  std::string rname(Process_Base::GenerateName(rampl));
  p_rproc=(*(*procs)[nlo_type::lo])[rname]->Get<Single_Process>();
  if (p_rproc==NULL) {
    msg_Debugging()<<"invalid real process '"<<rname<<"'\n";
    rampl->Delete();
    return Weight_Map();
  }
  Weight_Map ws;
  int stat(CalculateWeights(rampl,idmap,ws));
  rampl->Delete();
  if (stat==-1) ws.clear();
  return ws;
}

int CS_Gamma::CalculateWeights(Cluster_Amplitude *const ampl,
			       const std::map<size_t,size_t> &idmap,
			       Weight_Map &ws)
{
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *li(ampl->Leg(i));
    for (size_t j(Max((size_t)2,i+1));j<ampl->Legs().size();++j) {
      Cluster_Leg *lj(ampl->Leg(j));
      if (!p_rproc->Combinable(li->Id(),lj->Id())) continue;
      const Flavour_Vector &cf(p_rproc->CombinedFlavour(li->Id()|lj->Id()));
      for (size_t f(0);f<cf.size();++f) {
	for (size_t k(0);k<ampl->Legs().size();++k) {
	  Cluster_Leg *lk(ampl->Leg(k));
	  if (k==i || k==j) continue;
	  if (!CheckColors(li,lj,lk,cf[f])) continue;
	  if (!cf[f].Strong()) continue;/*QCD only so far*/
	  CS_Parameters cs(p_cluster->KT2(ampl,li,lj,lk,cf[f],p_ms));
	  if (cs.p_sf==NULL || !cs.p_sf->On()) continue;
	  Vec4D_Vector p(p_cluster->Combine(*ampl,i,j,k,cf[f],p_ms,cs.m_kin));
	  if (p.empty()) {
	    msg_Debugging()<<"combine failed for "<<ID(li->Id())<<"&"
			   <<ID(lj->Id())<<" <-> "<<ID(lk->Id())<<"\n";
	    continue;
	  }
	  ampl->SetKT2(cs.m_kt2);
	  ampl->SetZ(cs.m_z);
	  ampl->SetPhi(cs.m_phi);
	  ampl->SetKin(cs.m_kin);
	  Cluster_Amplitude *nampl(ampl->InitNext());
	  nampl->SetNIn(ampl->NIn());
	  nampl->SetMuF2(cs.m_kt2);
	  nampl->SetMuR2(ampl->MuR2());
	  nampl->Decays()=ampl->Decays();
	  nampl->SetProcs(ampl->Procs<void>());
	  nampl->SetDInfo(ampl->DInfo<void>());
	  Cluster_Leg *lijt(NULL), *lkt(NULL);
	  for (size_t l(0), m(0);l<ampl->Legs().size();++l) {
	    if (l==j) continue;
	    else if (l==i) {
	      nampl->CreateLeg(p[m],cf[f],ColorID(0,0),li->Id()|lj->Id());
	      nampl->Legs().back()->SetK(lk->Id());
	      lijt=nampl->Legs().back();
	    }
	    else {
	      Cluster_Leg *cl(ampl->Leg(l));
	      nampl->CreateLeg(p[m],cl->Flav(),cl->Col(),cl->Id());
	      if (cl==lk) lkt=nampl->Legs().back();
	    }
	    ++m;
	  }
	  if (CheckCore(nampl)) {
	    std::string pname(Process_Base::GenerateName(ampl));
	    const DInfo_Set &dinfo((*nampl->DInfo<DInfo_Map>())[pname]);
	    if (dinfo.find(DDip_ID(i,j,k))!=dinfo.end()) {
	      int stat(SingleWeight(nampl,li,lj,lk,cs,idmap,ws));
	      if (stat==-1) return -1;
	    }
	  }
	  ampl->DeleteNext();
	}
      }
    }
  }
  return 1;
}

int CS_Gamma::SingleWeight
(Cluster_Amplitude *const ampl,Cluster_Leg *const li,
 Cluster_Leg *const lj,Cluster_Leg *const lk,const CS_Parameters &cs,
 const std::map<size_t,size_t> &idmap,Weight_Map &ws)
{
#ifdef DEBUG__Trial_Weight
  DEBUG_FUNC(ID(li->Id())<<","<<ID(lj->Id())<<"<->"<<ID(lk->Id()));
#endif
  Splitting_Function_Base *cdip(cs.p_sf);
  if (cdip==NULL) return 0;
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"B config -> "<<*ampl<<" -> "<<cs<<" ( "
		 <<cdip->GetFlavourA()<<" -> "<<cdip->GetFlavourB()
		 <<" "<<cdip->GetFlavourC()<<" )\n";
  Cluster_Amplitude *pampl(ampl);
  while ((pampl=pampl->Next())!=NULL) msg_Debugging()<<*pampl<<"\n";
#endif
  cdip->SetFlavourSpec((lk->Id()&((1<<ampl->NIn())-1))?
		       lk->Flav().Bar():lk->Flav());
  double eta=1.0;
  if (cs.m_mode==1) eta=p_cluster->GetX(li,cdip)*cs.m_z;
  else if (cs.m_mode==2) eta=p_cluster->GetX(lk,cdip)*(1.0-cs.m_y);
  else if (cs.m_mode==3) eta=p_cluster->GetX(li,cdip)*cs.m_z;
  Weight_Value meps(Differential(ampl,0));
  meps.p_sf=cdip;
  meps.m_me*=cdip->SymFac();
#ifdef DEBUG__Trial_Weight
  double me=meps.m_me;
#endif
  meps.m_me*=(*cdip)(cs.m_z,cs.m_y,eta,cs.m_kt2,cs.m_q2,ampl)*
    cdip->MEPSWeight(cs.m_z,cs.m_y,eta,cs.m_kt2,cs.m_q2,ampl);
  meps.m_qij2=PDF::Qij2(li->Mom(),lj->Mom(),lk->Mom(),
			li->Flav(),lj->Flav());
  if (meps.m_me==0.0) {
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"zero matrix element\n";
#endif
    return 0;
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"add ( z = "<<cs.m_z<<", y = "<<cs. m_y
		 <<", kt = "<<sqrt(cs.m_kt2)<<" ) {\n  "<<*li
		 <<"\n  "<<*lj<<"\n  "<<*lk<<"\n} -> w = "
		 <<me<<" * "<<meps.m_me/me<<" -> "
		 <<meps.m_me<<" ( O = "<<cdip->RBMax()<<", S = "
		 <<cdip->SymFac()<<" )\n";
#endif
  ws[Weight_Key(idmap.find(li->Id())->second|idmap.find(lj->Id())->second,
		idmap.find(lk->Id())->second,cdip->GetFlavourA(),
		cdip->GetFlavourB(),cdip->GetFlavourC())]=meps;
  return 1;
}

bool CS_Gamma::CheckCore(const ATOOLS::Cluster_Amplitude *ampl) const
{
  std::vector<std::vector<size_t> > allids(ampl->Legs().size()-2);
  std::vector<std::set<std::pair<size_t,size_t> > > allncs(allids.size()-1);
  std::vector<size_t> &ids(allids.front());
  ids.resize(ampl->Legs().size());
  for (size_t i(0);i<ids.size();++i) ids[i]=ampl->Leg(i)->Id();
  for (size_t n(0);n<allids.size();) {
    std::set<std::pair<size_t,size_t> > &ncs(allncs[n]);
    std::vector<size_t> &ids(allids[n+1]=allids[n]);
    bool cb(false);
    for (std::vector<size_t>::iterator i(ids.begin());i!=ids.end();++i) {
      for (std::vector<size_t>::iterator j(i+1);j!=ids.end();++j) {
	if (p_rproc->Combinable(*i,*j) && 
	    ncs.find(std::pair<size_t,size_t>(*i,*j))==ncs.end()) {
	  ncs.insert(std::pair<size_t,size_t>(*i,*j));
	  *i|=*j;
	  j=ids.erase(j);
	  cb=true;
	  break;
	}
      }
      if (cb) break;
    }
    if (ids.size()==3) {
      if (p_rproc->Combinable(ids[0],ids[1]) &&
	  p_rproc->Combinable(ids[0],ids[2]) &&
	  p_rproc->Combinable(ids[1],ids[2])) {
#ifdef DEBUG__Trial_Weight
	msg_Debugging()<<"accept "<<Process_Base::GenerateName(ampl)<<"\n";
#endif
	return true;
      }
    }
    if (cb && ids.size()>3) ++n;
    else {
      if (n==0) break;
      ncs.clear();
      --n;
    }
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"reject "<<Process_Base::GenerateName(ampl)<<"\n";
#endif
  return false;
}

bool CS_Gamma::Reject()
{
  if (m_on==0) return false;
  Cluster_Amplitude *rampl=p_css->GetRealEmissionAmplitude();
  double wgt(TrialWeight(rampl));
  rampl->Delete();
  if (wgt>ran.Get()) {
    msg_Debugging()<<"w = "<<wgt<<" -> accept\n";
    return false;
  }
  msg_Debugging()<<"w = "<<wgt<<" -> reject\n";
  return true;
}

double CS_Gamma::TrialWeight(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  p_ms=ampl->MS();
  p_shower->SetMS(p_ms);
  Weight_Map ws(CalculateWeight(ampl));
  if (ws.empty()) return 1.0/p_shower->GetLast()[3]->SF()->RBMax();
  Parton *const *cur(p_shower->GetLast());
  double wgt(0.0);
  Weight_Value wact;
  Weight_Map::const_iterator ait;
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"Accumulate weights {\n";
#endif
  for (Weight_Map::const_iterator
	 wit(ws.begin());wit!=ws.end();++wit) {
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"  "<<wit->first<<" -> "<<wit->second;
#endif
    wgt+=wit->second.m_me;
    if ((wit->first.m_ij==(cur[0]->Id()|ampl->IdNew())) &&
	(wit->first.m_k==cur[2]->Id())) {
      ait=wit;
      wact=ait->second;
#ifdef DEBUG__Trial_Weight
      msg_Debugging()<<" <- active";
#endif
    }
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"\n";
#endif
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"} -> w = "<<wgt<<"\n";
#endif
  if (wact.m_me==-1.0)
    THROW(fatal_error,"No active splitting weight");
  ampl->SetMuF2(wact.m_muf2);
  ampl->SetMuR2(wact.m_mur2);
  double rme(wact.m_me/wgt*Differential(ampl,0).m_me);
  msg_Debugging()<<"me / ecss = "<<rme<<" / "<<wact.m_me
		 <<" = "<<rme/wact.m_me<<" ( O = "
		 <<wact.p_sf->RBMax()<<" )";
  RB_Data *rbd(wact.p_proc->Integrator()->RBMap()[ait->first.m_rbkey]);
  if (rbd->m_ktres>0.0) {
    ZH_Pair zh(p_css->ZHSplit(wact.m_b,wact.m_qij2,rbd));
    rme*=zh.first/(zh.first+zh.second);
    msg_Debugging()<<" -> "<<rme/wact.m_me;
  }
  msg_Debugging()<<"\n";
  return rme/(wact.m_me*wact.p_sf->RBMax());
}

void CS_Gamma::AddRBPoint(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC("");
  p_ms=ampl->MS();
  p_shower->SetMS(p_ms);
  Weight_Map ws(CalculateWeight(ampl));
  if (ws.empty()) return;
  double wsum(0.0);
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"Accumulate weights {\n";
#endif
  for (Weight_Map::const_iterator
	 wit(ws.begin());wit!=ws.end();++wit) {
#ifdef DEBUG__Trial_Weight
    msg_Debugging()<<"  "<<wit->first<<" -> "<<wit->second<<"\n";
#endif
    wsum+=wit->second.m_me;
  }
#ifdef DEBUG__Trial_Weight
  msg_Debugging()<<"} -> w = "<<wsum<<"\n";
#endif
  std::map<double,double> rmes;
  for (Weight_Map::const_iterator
	 bit(ws.begin());bit!=ws.end();++bit) {
    Process_Base *bproc(bit->second.p_proc);
    msg_Debugging()<<"Set weight for '"<<bproc->Name()
		   <<"'"<<bit->first<<" -> ";
    double wgt(bit->second.m_me);
    ampl->SetMuF2(bit->second.m_muf2);
    ampl->SetMuR2(bit->second.m_mur2);
    std::map<double,double>::iterator rit(rmes.find(ampl->MuF2()));
    double rme(rit!=rmes.end()?rit->second:Differential(ampl,0).m_me);
    rmes[ampl->MuF2()]=rme;
    rme*=wgt/wsum;
    msg_Debugging()<<"me / ecss = "<<rme<<" / "<<wgt
		   <<" = "<<rme/wgt;
    RB_Data *rbd(bproc->Integrator()->RBMap()[bit->first.m_rbkey]);
    if (rme/wgt>m_zhth && rbd->m_ktres==0.0) {
      double bmax(rbd->m_bmax);
      (*rbd)=RB_Data();
      rbd->m_bmax=bmax;
      rbd->m_ktres=m_ktres;
    }
    if (rbd->m_ktres>0.0) {
      // if ZH splitting get some statistics on m_bmax first
      if (rbd->m_n<2000) {
        rbd->AddPoint(1.0,0.0,bit->second.m_b);
        msg_Debugging()<<" -> not considering point, only taking stats on m_bmax.\n";
        return;
      }
      ZH_Pair zh(p_css->ZHSplit(bit->second.m_b,bit->second.m_qij2,rbd));
      rme*=zh.first/(zh.first+zh.second);
      msg_Debugging()<<" -> "<<rme/wgt;
    }
    msg_Debugging()<<"\n";
    rbd->AddPoint(rme/wgt,rme,bit->second.m_b);
    if (rme/wgt>1.0e2) {
      msg_Error()<<METHOD<<"(): (R/B)_ME/PS = "<<rme/wgt
                 <<" in the event #"<<rbd->m_n
		 <<" !\n  Event generation will hardly work for '"
		 <<bproc->Name()<<"'."<<std::endl;
    }
  }
}

Weight_Value CS_Gamma::Differential
(Cluster_Amplitude *const ampl,const int mode) const
{
#ifndef DEBUG__Differential
  int olv(msg->Level());
  msg->SetLevel(2);
#endif
  ProcessMap_Map *procs(ampl->Procs<ProcessMap_Map>());
  Process_Base::SortFlavours(ampl);
  std::string pname(Process_Base::GenerateName(ampl));
  Process_Map::const_iterator pit((*(*procs)[nlo_type::lo]).find(pname));
  if (pit==(*(*procs)[nlo_type::lo]).end()) 
    THROW(fatal_error,"Process '"+pname+"' not found");
  Weight_Value meps(pit->second);
  meps.m_b=meps.m_me=pit->second->Differential(*ampl,(mode&~1024)|1|2|4);
  meps.m_me*=pit->second->SymFac();
  if (mode&1024) {
    Process_Map::const_iterator sit((*(*procs)[nlo_type::rsub]).find(pname));
    if (sit==(*(*procs)[nlo_type::lo]).end()) 
      THROW(fatal_error,"Process '"+pname+"' not found");
    if (!pit->second->Selector()->JetTrigger
	(pit->second->Integrator()->Momenta(),
	 sit->second->GetSubevtList())) meps.m_me=0.0;
  }
  meps.m_muf2=ampl->MuF2();
  meps.m_mur2=ampl->MuR2();
#ifndef DEBUG__Differential
  msg->SetLevel(olv);
#endif
  return meps;
}

bool CS_Gamma::CheckColors
(const ATOOLS::Cluster_Leg *li,const ATOOLS::Cluster_Leg *lj,
 const ATOOLS::Cluster_Leg *lk,const ATOOLS::Flavour &mo) const
{
  if (mo.StrongCharge()==8) {
    if (lk->Flav().Strong()) return true;
  }
  else if (mo.Strong()) {
    if (lk->Flav().StrongCharge()==8 ||
	lk->Flav().StrongCharge()==-mo.StrongCharge()) return true;
  }
  else {
    if (lk->Flav().StrongCharge()!=8) return true;
  }
  return false;
}

namespace POWHEG {

  std::ostream &operator<<(std::ostream &str,const Weight_Key &k)
  {
    return str<<"["<<ATOOLS::ID(k.m_ij)<<","<<ATOOLS::ID(k.m_k)<<"]";
  }

  std::ostream &operator<<(std::ostream &str,const Weight_Value &w)
  {
    return str<<w.m_me<<"  "<<w.p_proc->Name()<<" [ "
	      <<w.p_sf->GetFlavourA()<<" -> "<<w.p_sf->GetFlavourB()
	      <<" "<<w.p_sf->GetFlavourC()<<" ] ( \\mu_F = "
	      <<sqrt(w.m_muf2)<<", \\mu_R = "<<sqrt(w.m_mur2)<<" ) ";
  }

}
