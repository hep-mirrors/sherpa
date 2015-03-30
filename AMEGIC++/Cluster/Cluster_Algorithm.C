#include "AMEGIC++/Cluster/Cluster_Algorithm.H"

#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "EXTRA_XS/Main/ME2_Base.H"
#include "AMEGIC++/Main/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>

using namespace AMEGIC;
using namespace PHASIC;
using namespace EXTRAXS;
using namespace ATOOLS;

Cluster_Algorithm::Cluster_Algorithm(ATOOLS::Mass_Selector *const ms):
  p_ms(ms), p_ampl(NULL), p_clus(NULL), p_combi(NULL) {}

Cluster_Algorithm::~Cluster_Algorithm()
{
  for (Flav_ME_Map::const_iterator xsit(m_xsmap.begin());
       xsit!=m_xsmap.end();++xsit) delete xsit->second;
  if (p_combi) delete p_combi;
}

bool Cluster_Algorithm::Cluster
(Process_Base *const xs,const size_t mode)
{
  DEBUG_FUNC("");
  p_proc=xs->GetReal();
  p_ampl=NULL;
  int nampl=p_proc->NumberOfDiagrams();
  int nlegs=p_proc->NIn()+p_proc->NOut();
  if (nampl==0) p_ct=NULL;
  else {
  Leg **legs(CreateLegs(nampl,nlegs));
  CreateTables(legs,nampl,mode);
  }
  ++m_cnt;
  if (nampl==0 || p_ct==NULL || p_ct->RScale()>0.0) {
    PHASIC::Process_Base *pb(p_proc->IsMapped()?
			     p_proc->MapProc():p_proc);
    double rscale((pb->Integrator()->Momenta()[0]+
		   pb->Integrator()->Momenta()[1]).Abs2());
    if (p_ct) rscale=p_ct->RScale();
    msg_Debugging()<<METHOD<<"(): {\n";
    p_ampl = Cluster_Amplitude::New();
    p_ampl->SetMS(p_ms);
    p_ampl->SetJF(p_proc->Selector()->GetSelector("Jetfinder"));
    p_ampl->SetNIn(p_proc->NIn());
    p_ampl->SetOrderEW(p_proc->MaxOrder(1));
    p_ampl->SetOrderQCD(p_proc->MaxOrder(0));
    std::vector<size_t> tids, atids;
    for (int i(0);i<pb->NIn()+pb->NOut();++i) {
      Flavour flav(i<pb->NIn()?p_proc->Flavours()[i].Bar():
		   p_proc->Flavours()[i]);
      Vec4D mom(i<pb->NIn()?-pb->Integrator()->Momenta()[i]:
		pb->Integrator()->Momenta()[i]);
      p_ampl->CreateLeg(mom,flav,ColorID(),1<<i);
      int sc(p_ampl->Legs().back()->Flav().StrongCharge());
      if (sc==0) p_ampl->Legs().back()->SetCol(ColorID(0,0));
      if (sc==3 || sc==8) {
	p_ampl->Legs().back()->SetCol(ColorID(Flow::Counter(),0));
	tids.push_back(i);
      }
      if (sc==-3 || sc==8) {
	p_ampl->Legs().back()->SetCol
	  (ColorID(sc==8?p_ampl->Legs().back()->Col().m_i:0,-1));
	atids.push_back(i);
      }
    }
    while (true) {
      std::random_shuffle(tids.begin(),tids.end(),*ran);
      size_t i(0);
      for (;i<tids.size();++i) if (tids[i]==atids[i]) break;
      if (i==tids.size()) break;
    }
    for (size_t i(0);i<atids.size();++i)
      p_ampl->Leg(atids[i])->SetCol
	(ColorID(p_ampl->Leg(atids[i])->Col().m_i,
		 p_ampl->Leg(tids[i])->Col().m_i));
    p_ampl->SetMuF2(pb->ScaleSetter()->Scale(stp::fac));
    p_ampl->SetMuR2(pb->ScaleSetter()->Scale(stp::ren));
    p_ampl->SetMuQ2(pb->ScaleSetter()->Scale(stp::res));
    PDF::CParam scale((p_proc->IsMapped()?p_proc->MapProc():p_proc)
                      ->ScaleSetter()->CoreScale(p_ampl));
    p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
    SetNMax(p_ampl,(1<<(p_proc->NIn()+p_proc->NOut()))-1,
            p_proc->Info().m_fi.NMaxExternal());

    // Ad-Hoc clustering for loop-induced gg -> non-coloured particles
    size_t np=p_proc->Flavours().size();
    bool non_col(true);
    for(size_t i(2); i<np-1; i++){
      if(p_proc->Flavours()[i].StrongCharge()){
	non_col = false;
	break;
      }
    }
    bool gg((p_proc->Flavours()[0].IsGluon() && p_proc->Flavours()[1].IsGluon() &&
	     p_proc->Flavours().back().IsGluon()) ||
	    (p_proc->Flavours()[0].IsGluon() && p_proc->Flavours()[1].IsQuark() &&
	     p_proc->Flavours().back().IsQuark()) ||
	    (p_proc->Flavours()[0].IsQuark() && p_proc->Flavours()[1].IsGluon() &&
	     p_proc->Flavours().back().IsQuark()));
    bool loop_ind((p_proc->MaxOrder(0)+p_proc->MaxOrder(1))==np);
    if (non_col && gg && loop_ind){
      ClusterSpecial4lLoop2();
    }
    // End Ad-Hoc clustering
    
    msg_Debugging()<<*p_ampl<<"\n";
    while (p_ampl->Prev()) {
      p_ampl=p_ampl->Prev();
      msg_Debugging()<<*p_ampl<<"\n";
    }
    msg_Debugging()<<"}\n";
    return true;
  }
  Convert();
  return true;
}

bool Cluster_Algorithm::FillLegs(Leg * alegs, Point * root, int & l, int maxl) 
{
  if (l>= maxl) {
    msg_Error()<<" Error in FillLegs() !!! "<<std::endl;
    return 0;
  }
  if (l==0) {
    size_t id(1<<root->number);
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    alegs[root->number].SetID(id);    
    l++;
  }
  if (root->left) {
    if (root->middle) return 0; // four vertex 
    return FillLegs(alegs,root->left,l,maxl)*FillLegs(alegs,root->right,l,maxl);
  } 
  else {
    size_t id(1<<root->number);
    alegs[root->number]=Leg(root);
    alegs[root->number].SetExternal(1);
    alegs[root->number].SetID(id);    
    l++;
    return 1;
  }
}

Leg **Cluster_Algorithm::CreateLegs(int &nampl,const int nlegs)
{
  Leg **legs(NULL);
  if (p_combi) delete p_combi;
  p_combi = 0;
  legs = new Leg *[nampl];
  for (int k=0;k<nampl;) {
    legs[k] = new Leg[nlegs];
    int l   = 0;
    if (FillLegs(legs[k],p_proc->Diagram(k),l,nlegs)) ++k;
    else {
      delete [] legs[k];
      --nampl;
    }
  }
  for (int k=0;k<nampl;++k) {
    for (int i(0);i<nlegs;++i) {
      Flavour fl(p_proc->Flavours()[i]);
      legs[k][i].SetMapFlavour(fl);
//       msg_Debugging()<<"set mapfl: "<<k<<", "<<i<<": "<<fl<<"\n";
    }
  }
  return legs;
}

void Cluster_Algorithm::CreateTables
(Leg ** legs,const int nampl,const size_t mode) 
{
  p_ct = 0;
  // if no combination table exist, create it
  int nin(p_proc->NIn()), nout(p_proc->NOut()), nlegs(nin+nout);
  Vec4D * amoms = new Vec4D[nlegs];
  for (int i=0;i<nin+nout;++i)  
    amoms[i]     = p_proc->Integrator()->Momenta()[i];
  if (!p_combi) {
    /*
      - copy moms to insert into Combine_Table (will be delete there)
      - create new Combine_Table with given momenta and given Jet-measure
      - initialise Combine_Table
      - determine best combination sheme
    */ 
    m_decids=p_proc->DecayInfos();
    p_combi = new Combine_Table(p_proc,p_ms,p_clus,amoms,0,&m_decids);
    p_combi->FillTable(legs,nlegs,nampl);   
    p_ct = p_combi->CalcJet(nlegs,NULL,mode,(mode&512)?1:((mode&16384)?1:0)); 
    if (p_ct==NULL && !(mode&512)) {
      msg_Debugging()<<"trying unordered configuration (top level)\n";
      p_ct = p_combi->CalcJet(nlegs,NULL,mode,0); 
    }
  }
  else {
    // use the existing combination table and determine best combination sheme
    p_ct = p_combi->CalcJet(nlegs,amoms,mode,(mode&512)?1:((mode&16384)?1:0));
    if (p_ct==NULL && !(mode&512)) {
      msg_Debugging()<<"trying unordered configuration (top level)\n";
      p_ct = p_combi->CalcJet(nlegs,NULL,mode,0); 
    }
  }
  //  delete [] amoms;
}

void Cluster_Algorithm::Convert()
{
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
  Selector_Base *jf=p_proc->Selector()
    ->GetSelector("Jetfinder");
  Combine_Table *ct_tmp(p_ct);
  while (ct_tmp->Up()) ct_tmp=ct_tmp->Up();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(jf);
  p_ampl->SetNIn(p_proc->NIn());
  p_ampl->SetOrderEW(p_proc->MaxOrder(1));
  p_ampl->SetOrderQCD(p_proc->MaxOrder(0));
  PHASIC::Process_Base *pb(p_proc->IsMapped()?
			   p_proc->MapProc():p_proc);
  double muf2(pb->ScaleSetter()->Scale(stp::fac));
  double mur2(pb->ScaleSetter()->Scale(stp::ren));
  double muq2(pb->ScaleSetter()->Scale(stp::res));
  for (int i(0);i<ct_tmp->NLegs();++i) {
    size_t id(ct_tmp->GetLeg(i).ID());
    Flavour flav(i<pb->NIn()?ct_tmp->Flav(i).Bar():ct_tmp->Flav(i));
    Vec4D mom(i<pb->NIn()?-ct_tmp->Momentum(i):ct_tmp->Momentum(i));
    p_ampl->CreateLeg(mom,flav,ColorID(0,0),id);
  }
  p_ampl->SetMuQ2(muq2);
  p_ampl->SetMuR2(mur2);
  p_ampl->SetMuF2(muf2);
  Cluster_Amplitude *eampl(p_ampl);
  while (ct_tmp->Down()) {
    int iwin, jwin, kwin, kmode;
    double mu2;
    double kt2qcd(ct_tmp->GetWinner(iwin,jwin,kwin,mu2,kmode));
    if (iwin>jwin) std::swap<int>(iwin,jwin);
    ct_tmp=ct_tmp->Down();
    const Leg &win(ct_tmp->GetLeg(iwin));
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl=p_ampl->InitNext();
    p_ampl->SetMS(p_ms);
    p_ampl->SetJF(jf);
    p_ampl->SetNIn(ampl->NIn());
    ampl->SetKT2(kt2qcd);
    ampl->SetMu2(mu2);
    for (int i(0);i<ct_tmp->NLegs();++i) {
      size_t id(ampl->Leg(i<jwin?i:i+1)->Id());
      Flavour flav(i<pb->NIn()?ct_tmp->Flav(i).Bar():ct_tmp->Flav(i));
      Vec4D mom(i<pb->NIn()?-ct_tmp->Momentum(i):ct_tmp->Momentum(i));
      if (i==iwin) id+=ampl->Leg(jwin)->Id();
      p_ampl->CreateLeg(mom,flav,ColorID(0,0),id);
      p_ampl->Legs().back()->SetStat(1);
      if (i==iwin) {
	p_ampl->Legs().back()->SetK(ampl->Leg(kwin)->Id());
	ampl->SetIdNew(ct_tmp->Up()->GetLeg(jwin).ID());
	if (win.Point()->t>10) {
	  size_t dmax(win.Point()->t>10?win.Point()->t-10:0);
	  if (dmax==0) dmax=IdCount(id);
	  p_ampl->Legs().back()->SetStat
	    (p_ampl->Legs().back()->Stat()|2);
	  SetNMax(p_ampl->Prev(),id,dmax);
	}
	if (kmode)
	  p_ampl->Legs().back()->SetStat
	    (p_ampl->Legs().back()->Stat()|4);
      }
    }
    p_ampl->SetMuQ2(ampl->MuQ2());
    p_ampl->SetMuR2(ampl->MuR2());
    p_ampl->SetMuF2(ampl->MuF2());
    p_ampl->Decays()=ct_tmp->Decays();
    p_ampl->SetOrderEW(ampl->OrderEW()-win.OrderQED());
    p_ampl->SetOrderQCD(ampl->OrderQCD()-win.OrderQCD());
    p_ampl->SetKin(win.Kin());
  }
  p_ampl->SetProc(p_proc);
  PDF::CParam scale((p_proc->IsMapped()?p_proc->MapProc():p_proc)
		    ->ScaleSetter()->CoreScale(p_ampl));
  p_ampl->SetKT2(scale.m_kt2);
  p_ampl->SetMu2(scale.m_mu2);
  size_t nmax(p_proc->Info().m_fi.NMaxExternal());
  p_ampl->Decays()=p_proc->Info().m_fi.GetDecayInfos();
  SetNMax(p_ampl,(1<<(p_proc->NIn()+p_proc->NOut()))-1,nmax);
  while (p_ampl->Prev()) {
    msg_Debugging()<<*p_ampl<<"\n";
    p_ampl=p_ampl->Prev();
  }
  msg_Debugging()<<*p_ampl<<"\n";
  msg_Debugging()<<"}\n";
}

void Cluster_Algorithm::SetNMax(Cluster_Amplitude *const ampl,
				const size_t &id,const size_t &nmax) const
{
  if (ampl==NULL) return;
  for (size_t i(0);i<ampl->Legs().size();++i) {
    Cluster_Leg *cli(ampl->Leg(i));
    if (cli->Id()&id) {
      cli->SetNMax(nmax);
      if (cli->Stat()!=3) 
	SetNMax(ampl->Prev(),cli->Id(),nmax);
    }
  }
}

void Cluster_Algorithm::ClusterSpecial4lLoop2()
{
  DEBUG_FUNC(*p_ampl);
  size_t emitted_idx=p_proc->Flavours().size()-1;
      
  PDF::CParam c0=p_clus->KPerp2(*p_ampl, 0, emitted_idx, 1, Flavour(kf_gluon), p_ms);
  PDF::CParam c1=p_clus->KPerp2(*p_ampl, 1, emitted_idx, 0, Flavour(kf_gluon), p_ms);

  int winner=0;
  PDF::CParam win=c0;
  if (ran->Get()*(1.0/c0.m_op2+1.0/c1.m_op2)>1.0/c0.m_op2) {
    winner=1;
    win=c1;
  }
      
  Cluster_Amplitude *ampl(p_ampl);
  p_ampl=p_ampl->InitNext();
  p_ampl->SetMS(p_ms);
  p_ampl->SetJF(p_proc->Selector()->GetSelector("Jetfinder"));
  p_ampl->SetNIn(p_proc->NIn());
  p_ampl->SetOrderEW(p_proc->MaxOrder(1));
  p_ampl->SetOrderQCD(p_proc->MaxOrder(0)-1);
  ampl->SetKT2(win.m_kt2);
  ampl->SetMu2(win.m_mu2);

  Vec4D_Vector clustered_moms=p_clus->Combine
    (*ampl, winner, emitted_idx, 1-winner, Flavour(kf_gluon),p_ms);

  int color[2] = { Flow::Counter(), Flow::Counter() };
  for (int i=0;i<2; ++i) {
    size_t id=(1<<i);
    if (i==winner) id+=(1<<emitted_idx);
    p_ampl->CreateLeg(clustered_moms[i],Flavour(kf_gluon),ColorID(color[i],color[1-i]),id);
    p_ampl->Legs().back()->SetStat(1);
    p_ampl->Legs().back()->SetNMax(p_proc->Info().m_fi.NMaxExternal());
    if (i==winner) {
      p_ampl->Legs().back()->SetK(1<<(1-i));
    }
  }
  for (int i=2;i<p_proc->Flavours().size()-1;++i) {
    p_ampl->CreateLeg(clustered_moms[i],p_proc->Flavours()[i],ColorID(0,0),(1<<i));
    p_ampl->Legs().back()->SetNMax(p_proc->Info().m_fi.NMaxExternal());
    p_ampl->Legs().back()->SetStat(1);
  }
  p_ampl->SetKin(win.m_kin);
  p_ampl->SetMuR2(ampl->MuR2());
  p_ampl->SetMuF2(ampl->MuF2());
  p_ampl->SetMuQ2(ampl->MuQ2());
  (p_proc->IsMapped()?p_proc->MapProc():p_proc)->ScaleSetter()->CoreScale(p_ampl);
}
