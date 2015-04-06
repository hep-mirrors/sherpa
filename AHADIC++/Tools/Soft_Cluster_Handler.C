#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Return_Value.H"

using namespace AHADIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(bool ana) :
  p_singletransitions(hadpars->GetSingleTransitions()), 
  p_doubletransitions(hadpars->GetDoubleTransitions()),
  p_as((Strong_Coupling*)s_model->GetScalarFunction("strong_cpl")),
  m_transoffset(hadpars->Get(std::string("Offset_C->H"))), 
  m_decayoffset(hadpars->Get(std::string("Offset_C->HH"))), 
  m_kappa(hadpars->Get(std::string("MassExponent_C->H"))), 
  m_lambda(hadpars->Get(std::string("WidthExponent_C->H"))), 
  m_chi(hadpars->Get(std::string("MassExponent_C->HH"))), 
  m_pt2max(sqr(hadpars->Get(string("ptmax")))),
  m_pt2maxfac(sqr(hadpars->Get(std::string("ptmax_factor")))),
  m_pt02(hadpars->Get(std::string("pt02"))), 
  m_transitions(0), m_dtransitions(0), m_decays(0), 
  m_forceddecays(0), m_lists(0), m_update(0),
  m_ana(ana)
{
  if (m_ana) {
    m_histograms[string("PT_HH")]  = new Histogram(0,0.,10.,100);
    m_histograms[string("PT2_HH")] = new Histogram(0,0.,100.,2000);
    m_histograms[string("MassTransition")]       = new Histogram(0,0.,8.,100);
    m_histograms[string("HadronMassTransition")] = new Histogram(0,0.,8.,100);
  }
}

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{
  msg_Tracking()<<"@@@ "<<METHOD<<": "
		<<m_transitions<<" transitions, "
		<<m_dtransitions<<" double transitions, "
		<<m_decays<<" decays,\n"
		<<m_update<<" calls to UpdateTransitions, and "
		<<m_forceddecays<<" forced decays.\n"
		<<"@@@ "<<METHOD<<": "
		<<m_lists<<" transitions from original dipole list.\n";
  if (m_ana) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}

bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin, Blob * blob)
{
  if (clin->size()==1) {
    //msg_Out()<<METHOD<<" for cluster list with one cluster only!\n";
    if (TreatSingleCluster(*clin->begin(),true))
      return AttachHadronsToBlob(clin,blob);
    if (EnforcedDecay(*clin->begin(),blob,true,clin,false)) {
      AttachHadronsToBlob(clin,blob);
      //msg_Out()<<METHOD<<" enforces decay of cluster:\n"<<(*blob)<<",\n"
	      // <<"   cluster list has "<<clin->size()<<" elements left.\n";
      return true;
    }
  }
  //msg_Out()<<METHOD<<" for "<<clin->size()<<" clusters.\n";
  if (!CheckListForTreatment(clin)) return true;
  m_lists++;
  double E(-1.);
  //msg_Out()<<METHOD<<" checking if all decays/transitions are allowed.\n";
  if (!CheckIfAllowed(clin,E)) {
    do {
      if (!UpdateTransitions(clin)) {
	msg_Error()<<"Error in "<<METHOD<<", E = "<<E<<" : \n"
		   <<"   Cannot update any further.  "
		   <<"Leave and hope for the best.\n"
		   <<(*clin)<<"\n";
	return false;
      }
    } while (!CheckIfAllowed(clin,E));
    if (!EnforcedTransition(clin)) {
      msg_Error()<<"Error in "<<METHOD<<", E = "<<E<<" : \n"
		 <<"   Cannot enforce any meaningful transition. "
		 <<"Leave and hope for the best.\n"
		 <<(*clin)<<"\n";
      return false;
    }
  }
  return (ShiftMomenta(clin) && AttachHadronsToBlob(clin,blob));
}

bool Soft_Cluster_Handler::TreatSingleCluster(Cluster * cluster,bool mustdecay)
{
  switch (CheckCluster(cluster,true,mustdecay)) {
  case 2:  break;
  case 1:  cluster->push_back(Flavour(kf_photon)); break;
  case 0:  break;
  case -1: return false;
  default: break;
  }
  //msg_Out()<<METHOD<<";\n"<<(*cluster)<<"\n";
  return true;
}

bool Soft_Cluster_Handler::TreatClusterDecay(Cluster_List * clin, Blob * blob)
{
  if (!CheckListForTreatment(clin)) {
    return true;
  }
  Cluster_Iterator cit(clin->begin());
  Cluster * left((*cit++)), * right((*cit));
  Cluster * cluster(right->GetPrev());
  msg_Debugging()<<METHOD
		 <<"------------------------------------------------------\n"
		 <<" for "<<clin->size()<<" clusters with masses "
		 <<left->Mass()<<" + "<<right->Mass()<<", "
		 <<"original cluster: mass = "<<cluster->Mass()<<".\n";

  double E(-1.);
  while (!CheckIfAllowed(clin,E)) {
    if (UpdateTransitions(clin)) continue;
    if (left->GetPrev()!=right->GetPrev()) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<clin->size()<<" clusters) : \n"
		 <<"   No common previous cluster.\n";
      return false;
    }
    clin->clear();
    //msg_Out()<<METHOD<<" --> EnforcedDecay for \n"<<(*cluster)<<"\n";
    return EnforcedDecay(cluster,blob,true,clin,false);
    // probably better:   
    // if (!EnforcedDecay(cluster,blob,true,clin) return false;
  }  
  return (UpdateClusterDecayKinematics((*clin->begin())->GetPrev()) &&
	  AttachHadronsToBlob(clin,blob));
}

bool Soft_Cluster_Handler::
AttachHadronsToBlob(Cluster_List * clin,Blob * blob)
{
  Cluster_Iterator cit(clin->begin());
  Particle * part;
  Cluster * cluster;
  while (cit!=clin->end()) {
    cluster = (*cit);
    //msg_Out()<<METHOD<<" for "<<cluster->size()<<".\n";
    switch (cluster->size()) {
    case 1:
      part = cluster->GetSelf();
      part->SetFinalMass();
      blob->AddToOutParticles(part);
      msg_Tracking()<<"$$ attach one hadron ("<<part->Flav()<<", "
		    <<part->Momentum()<<", "
		    <<"pt = "<<part->Momentum().PPerp()<<", "
		    <<"y = "<<part->Momentum().Y()<<") "
		    <<"from cluster "<<cluster->Number()<<", "
		    <<"m = "<<cluster->Mass()<<".\n";
      delete cluster->GetTrip();
      delete cluster->GetAnti();
      delete cluster;
      cit = clin->erase(cit);
      break;
    case 2:
      FixHHDecay(cluster,blob,(*cluster)[0],(*cluster)[1]);
      delete cluster->GetTrip();
      delete cluster->GetAnti();
      delete cluster;
      cit = clin->erase(cit);
      break;      
    case 0:
      cit++;
      break;
    default:
      cit++;
      break;
    }
  }
  return true;
}

bool Soft_Cluster_Handler::UpdateClusterDecayKinematics(Cluster * cluster)
{  
  size_t nclusters(cluster->GetClusters()->size());
  if (nclusters==0) return true;
  vector<Vec4D>  moms;
  vector<double> masses;
  bool stretchthem(false);
  for (Cluster_Iterator cit=cluster->GetClusters()->begin();
       cit!=cluster->GetClusters()->end();cit++) {
    if ((*cit)->size()==1) {
      stretchthem=true;
      break;
    }
  }
  if (stretchthem) {
    cluster->BoostInCMS();
    for (Cluster_Iterator cit=cluster->GetClusters()->begin();
	 cit!=cluster->GetClusters()->end();cit++) {
      moms.push_back((*cit)->Momentum());
      masses.push_back((*cit)->size()!=1?(*cit)->Mass():(**cit)[0].HadMass());
    }
    Momenta_Stretcher stretcher;
    if (!stretcher.ZeroThem(0,moms) ||
	!stretcher.MassThem(0,moms,masses)) return false;
    size_t i(0);
    for (Cluster_Iterator cit=cluster->GetClusters()->begin();
	 cit!=cluster->GetClusters()->end();cit++) {
      (*cit)->SetMomentum(moms[i]);
      if ((*cit)->size()!=1) (*cit)->RescaleMomentum(moms[i]);
      i++;
    }
    cluster->BoostBack();
  }
  return true;
}

bool Soft_Cluster_Handler::CheckListForTreatment(Cluster_List * clin) {
  Cluster_Iterator cit;
  Cluster * cluster;
  int    size(0),count(0);
  for (cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster==NULL || !cluster->Active()) continue;
    size = CheckCluster(cluster,false);
    count += size;
  }
  if (count==0) return false;
  return true;
}


int Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool lighter,
				       bool mustdecay)
{
  cluster->clear();

  Flavour haddec1(Flavour(kf_none)), haddec2(Flavour(kf_none));
  Flavour hadtrans(Flavour(kf_none));
  double decayweight(DecayWeight(cluster,haddec1,haddec2,mustdecay));
  double transformweight(TransformWeight(cluster,hadtrans,lighter,mustdecay));
  //msg_Out()<<METHOD<<"("<<mustdecay<<lighter<<") yields "
	  // <<"dec = "<<decayweight<<"("<<haddec1<<" "<<haddec2<<"), "
	  // <<"trans = "<<transformweight<<"("<<hadtrans<<").\n"
	  // <<(*cluster)<<"\n";
	   
  if (mustdecay) {
    if (decayweight>0.) {
      cluster->push_back(haddec1);
      cluster->push_back(haddec2);
      m_decays++;
      return 2;
    }
    else if (transformweight>0.) {
      cluster->push_back(hadtrans);
      cluster->push_back(Flavour(kf_photon));
      m_decays++;
      return 2;
    }
    cluster->clear();
    return -1;
  }
  if (decayweight>0. || transformweight>0.) {
    double totweight((Max(0.,decayweight)+Max(0.,transformweight))*0.9999999);
    if (decayweight>0. && 
	(totweight*ran->Get()>decayweight || transformweight==0.)) {
      m_decays      += 1;
      cluster->push_back(haddec1);
      cluster->push_back(haddec2);
      //msg_Out()<<METHOD<<" yields decay to "<<haddec1<<" + "<<haddec2<<".\n";
      return 2;
    }
    if (transformweight>0.) {
      m_transitions += 1;
      cluster->push_back(hadtrans);
      //msg_Out()<<METHOD<<" yields transition to "<<hadtrans<<".\n";
      return 1;
    }
  }
  else {
    if (transformweight<=0. && hadtrans!=Flavour(kf_none)) {
      transformweight = TransformWeight(cluster,hadtrans,false,true);
      //msg_Out()<<METHOD<<" tries to force it for\n"<<(*cluster)<<"\n";
      //msg_Out()<<"Gotcha: "<<transformweight<<", "<<hadtrans<<"!\n";
      m_transitions += 1;
      cluster->push_back(hadtrans);
      return 1;
    }
  }
  cluster->clear();
  return 0;
}

bool Soft_Cluster_Handler::CheckIfAllowed(Cluster_List * clin,double & E) {
  double totmass(0.);
  Vec4D  totmom(0.,0.,0.,0.);
  Cluster * cluster;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    msg_Tracking()<<METHOD<<"("<<cluster->Number()<<").\n";
    switch (cluster->size()) {
    case 1: 
      totmass += (*cluster)[0].HadMass();
      break;
    case 2:
    case 0:
    default:
      totmass += cluster->Mass();
      break;
    }
    if (E<0) totmom += cluster->Momentum();
  }
  if (E<0) E = sqrt(totmom.Abs2());
  return (totmass<E);
}

bool Soft_Cluster_Handler::UpdateTransitions(Cluster_List * clin) {
  m_update++;
  Cluster * cluster, * winner(NULL);
  Flavour hadron1,hadron2,winhad1,winhad2;
  double  wt, maxwt(0.);
  int     winno(0);
  bool    found(false);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    wt      = -1.;
    msg_Tracking()<<METHOD<<"("<<cluster->Number()<<").\n";
    if (cluster->size()==1) {
      hadron1 = (*cluster)[0];
      wt      = TransformWeight(cluster,hadron1,true,true);
      if (wt>maxwt && 
	  (cluster->size()==0 || 
	   (cluster->size()==1 && hadron1!=(*cluster)[0]))) {
	winner  = cluster;
	winhad1 = hadron1;
	winno   = 1;
	maxwt   = wt;
	found   = true; 
      }
    }
    else if (wt<0. || cluster->size()==2) {
      wt     = DecayWeight(cluster,hadron1,hadron2,true);
      if (wt>maxwt && 
	  (cluster->size()==0 || 
	   (cluster->size()==2 && 
	    hadron1!=(*cluster)[0] && hadron2!=(*cluster)[1]))) {
	winner  = cluster;
	winhad1 = hadron1;
	winhad2 = hadron2;
	winno   = 2;
	maxwt   = wt;
	found   = true; 
      }
    }
  }
  if (found) {
    if (winno==1) {
      winner->clear();
      winner->push_back(winhad1);
      return true;
    }
    else if (winno==2) {
      winner->clear();
      winner->push_back(winhad1);
      winner->push_back(winhad2);
      return true;
    }
  }
  return false;
}

bool Soft_Cluster_Handler::
ClusterAnnihilation(Cluster * cluster,Flavour & had1,Flavour & had2) {
  int kfc1(int(cluster->GetTrip()->m_flav.Kfcode())); 
  int kfc2(int(cluster->GetAnti()->m_flav.Kfcode())); 
  kf_code kfc11(kfc1/1000),kfc12((kfc1-kfc11*1000)/100);
  kf_code kfc21(kfc2/1000),kfc22((kfc2-kfc21*1000)/100);
  Flavour fl1(kfc11), fl2(kfc12), fl3(kfc21), fl4(kfc22);
  //msg_Out()<<METHOD<<"("<<fl1<<", "<<fl2<<") & ("<<fl3<<", "<<fl4<<").\n";
  fl1 = fl1.Bar();
  fl2 = fl2.Bar();
  Proto_Particle *pp1(new Proto_Particle(fl1,cluster->GetTrip()->m_mom/2.,'l'));
  Proto_Particle *pp2(new Proto_Particle(fl2,cluster->GetTrip()->m_mom/2.,'l'));
  Proto_Particle *pp3(new Proto_Particle(fl3,cluster->GetAnti()->m_mom/2.,'l'));
  Proto_Particle *pp4(new Proto_Particle(fl4,cluster->GetAnti()->m_mom/2.,'l'));
  bool order(ran->Get()>0.5?true:false);
  Cluster cluster1((order?pp3:pp4),pp1), cluster2((order?pp4:pp3),pp2);
  Flavour_Pair pair1, pair2;
  pair1.first = cluster1.GetTrip()->m_flav;
  pair1.first = cluster1.GetAnti()->m_flav;
  pair2.first = cluster2.GetTrip()->m_flav;
  pair2.first = cluster2.GetAnti()->m_flav;
  double mass(cluster->Mass());
  double wt1(TransformWeight(&cluster1,had1,true,true));
  double wt2(TransformWeight(&cluster2,had2,true,true));
  if (wt1<0. || wt2<0.) return false;
  bool lighter1(wt1>0.), lighter2(wt2>0.);
  while (had1.Mass()+had2.Mass()>mass && lighter1 && lighter2) {
    lighter1 = wt1>0.;
    lighter2 = wt2>0.;
    if (wt1<=wt2) {
      if (lighter1)      wt1 = TransformWeight(&cluster1,had1,true,true);
      else if (lighter2) wt2 = TransformWeight(&cluster2,had2,true,true);
      else return false;
    }
    else {
      if (lighter2)      wt2 = TransformWeight(&cluster2,had2,true,true);
      else if (lighter1) wt1 = TransformWeight(&cluster1,had1,true,true);
      else return false;
    }
  }
  return true;
}

bool Soft_Cluster_Handler::
EnforcedDecay(Cluster * cluster, Blob * blob,const bool & constrained,
	      Cluster_List * clin,const bool & attach) {
  Flavour had1,had2;
  double weight(DecayWeight(cluster,had1,had2,true)), weight1(-1.);
  if (weight<=0.) {
    if (cluster->GetTrip()->m_flav.IsDiQuark() && 
	cluster->GetAnti()->m_flav.IsDiQuark()) {
      //msg_Out()<<METHOD<<": must annihilate cluster.\n";
      if (!ClusterAnnihilation(cluster,had1,had2)) return false;
    }
    else {
      weight1 = TransformWeight(cluster,had1,true,true);
      if (weight1<=0.) {
	msg_Tracking()<<"Error in "<<METHOD<<" : \n"
		      <<"   No suitable single transition found, "
		      <<"will return false and hope for the best.\n";
	return false;
      }
      had2 = Flavour(kf_photon);
    }
    m_forceddecays++;m_decays--;
  }

  //msg_Out()<<METHOD<<" --> "<<had1<<" + "<<had2<<"\n";
  cluster->push_back(had1);
  cluster->push_back(had2);
  if (attach) {
    FixHHDecay(cluster,blob,had1,had2,constrained);
  }
  return true;
}

bool Soft_Cluster_Handler::EnforcedTransition(Cluster_List * clin) {
#ifdef AHAmomcheck
  Vec4D checkbef(SumMomentum(clin));
#endif
  size_t size(0);
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  Cluster * cluster;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    switch (cluster->size()) {
    case 1: 
      masses.push_back((*cluster)[0].HadMass());
      momenta.push_back(cluster->Momentum());
      ++size;
      break;
    case 2:
    case 0:
      masses.push_back(sqrt(Max(cluster->GetTrip()->m_mom.Abs2(),0.)));
      masses.push_back(sqrt(Max(cluster->GetAnti()->m_mom.Abs2(),0.)));
      momenta.push_back(cluster->GetTrip()->m_mom);
      momenta.push_back(cluster->GetAnti()->m_mom);
      size+=2;
    default:
      break;
    }
  }
  if (!hadpars->AdjustMomenta(size,&momenta.front(),&masses.front())) {
    if (size>1 /*&& msg->LevelIsDebugging()*/) {
      msg_Tracking()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : \n"
		    <<"   Could not adjust momenta for : \n";
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Tracking()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Tracking()<<" ("<<((**cit)[0])<<" )";
	msg_Tracking()<<" for \n"<<(**cit)<<"\n";
	i++;
      }
      msg_Tracking()<<"   Will possibly lead to retrying the event.\n";
    }
    return false;
  }
  int pos(0);

#ifdef AHAmomcheck
  Vec4D checkaft(SumMomentum(clin));
#endif
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    switch (cluster->size()) {
    case 1:
      cluster->SetFlav((*cluster)[0]);
      cluster->SetMomentum(momenta[pos]);
      break;
    case 2: 
#ifdef AHAmomcheck
      checkaft += momenta[pos];
#endif
      cluster->GetTrip()->m_mom=momenta[pos++];
      cluster->GetAnti()->m_mom=momenta[pos];
      cluster->Update();		
      if (cluster->Mass()<(*cluster)[0].HadMass()+(*cluster)[1].HadMass()) {
	cluster->clear();
      }
      break;
    case 0:
      cluster->GetTrip()->m_mom=momenta[pos++];
      cluster->GetAnti()->m_mom=momenta[pos];
      cluster->Update();		
      break;
    default:
      break;
    }
    pos++;
  }
#ifdef AHAmomcheck
  double Q2(dabs((checkbef-checkaft).Abs2()));
  if (Q2>1.e-12 || IsNan(Q2)) {
    msg_tracking()<<METHOD<<" yields a momentum violation for  "<<size<<" : \n"
		  <<"   "<<checkbef<<" - "<<checkaft<<" --> "
		  <<(checkbef-checkaft).Abs2()<<"("<<size<<").\n"
		  <<(*clin)<<"\n";
  }
  else msg_Tracking()<<METHOD<<" satisfied four-momentum conservation.\n";
#endif
  m_transitions += 1;
  return true;  
}


double Soft_Cluster_Handler::
TransformWeight(Cluster * cluster,Flavour & hadron,
		const bool & lighter,const bool & enforce)
{
  Flavour_Pair fpair;
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;
  if (p_singletransitions->GetTransitions()->find(fpair)==
      p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Error in "<<METHOD<<" for cluster\n"
	       <<(*cluster)
	       <<"   illegal flavour combination.\n"
	       <<"   Will return 0 and hope for the best.\n";
    return 0.;
  }

  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;

  double MC(cluster->Mass());
  double critM(p_singletransitions->GetLightestMass(fpair)*(1.-m_transoffset)+
	       p_singletransitions->GetHeaviestMass(fpair)*m_transoffset);
  //msg_Out()<<METHOD<<" for M = "<<MC<<" (crit = "<<critM<<")\n";
  if (!enforce && MC>critM) {
    hadron = Flavour(kf_none);
    return 0.;
  }
  Single_Transition_Miter stiter = 
    p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    hadron = Flavour(kf_none);
    return -1.;
  }
  Single_Transition_List * stl(stiter->second);
  Single_Transition_Siter  start(stl->begin()),siter;
  if (lighter) {
    if (hadron!=Flavour(kf_none)) {
        siter=stl->find(hadron);
        if (siter == stl->end()) {
            return 0.;
        }
        else {
            start=siter;
            start++; // Jump to the next lightest hadron
        }
    }
    else {
      for (siter=start;siter!=stl->end();siter++) {
	if (siter->first.Mass()<MC) {
	  start=siter;
	  break;
	}
      }
    }
  }

  double wt(0.),totweight(0.);
  
  for (siter=start;siter!=stl->end();siter++) {
    if (siter->first!=hadron) {
      wt  = TransformKin(MC,siter->first,enforce);
      wt *= siter->second;
    }
    totweight += wt;
  }
  //msg_Out()<<METHOD<<" with totweight = "<<totweight<<".\n";
//  for (siter=start;siter!=stl->end();siter++) {
 //   msg_Out()<<siter->first<<" yields "<<siter->second<<"  --> "
//	     <<TransformKin(MC,siter->first,enforce)<<".\n";
  //}

  double disc(totweight * 0.9999999999*ran->Get());
  for (siter=start;siter!=stl->end();siter++) {
    if (siter->first!=hadron) {
      wt  = TransformKin(MC,siter->first,enforce);
      wt *= siter->second;
    }
    disc -= wt;
    if (disc<=0.) {
      hadron = siter->first;
      disc   = wt;
      break;
    }
  }
  return wt/(16.*M_PI*MC);
}

double Soft_Cluster_Handler::
TransformKin(const double MC,const Flavour & flav,const bool & enforce) {
  double mass2(sqr(flav.HadMass()));
  double width2(sqr(Max(flav.Width(),1.e-6)));
  if (!enforce && sqr(MC*MC-mass2)>10.*mass2*width2) return 0.;
  return
    pow(sqr(mass2)/(sqr(MC*MC-mass2) + mass2*width2),m_kappa) * 
    pow(mass2*width2/(sqr(MC*MC-mass2) + mass2*width2),m_lambda);
}

double Soft_Cluster_Handler::
DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2,
	    const bool & enforce)
{
  Flavour_Pair flpair;
  flpair.first  = cluster->GetTrip()->m_flav;
  flpair.second = cluster->GetAnti()->m_flav;
  Double_Transition_Miter dtliter = 
    p_doubletransitions->GetTransitions()->find(flpair);
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg_Error()<<"Error in "<<METHOD<<" for cluster\n"
	       <<(*cluster)
	       <<"   illegal flavour combination.\n"
	       <<"   Will return 0 and hope for the best.\n";
    return 0.;
  }
  double MC(cluster->Mass()),MC2(MC*MC);
  
  if (p_doubletransitions->GetLightestMass(flpair)>MC && enforce) {
    return -1.;
  }

  bool hit1(cluster->GetTrip()->m_info=='L'||cluster->GetTrip()->m_info=='B');
  bool hit2(cluster->GetAnti()->m_info=='L'||cluster->GetAnti()->m_info=='B');
  bool direct(false);
  hit1 = hit2 = false;
  if (hit1 && hit2 &&
      (flpair.first==Flavour(kf_b) || flpair.first==Flavour(kf_c)) &&
      (flpair.first==Flavour(kf_b) || flpair.first==Flavour(kf_c))) {
    if (flpair.first.HadMass()+flpair.second.HadMass() >
	ran->Get()*(MC-flpair.first.HadMass()-flpair.second.HadMass()))
      direct = true;
  }
  double critM(p_doubletransitions->GetLightestMass(flpair)*(1.-m_decayoffset)+
	       p_doubletransitions->GetHeaviestMass(flpair)*m_decayoffset);
  if (!direct && !enforce && MC>critM) {
    had1 = had2 = Flavour(kf_none);
    return 0.;
  }
  double totweight(0.),m1,m2,wt(1.),wfweight(0.),wfmax(0.);
  double tm(cluster->GetTrip()->m_flav.HadMass());
  double am(cluster->GetTrip()->m_flav.HadMass());
  Flavour max1, max2;
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.HadMass();
    m2  = decit->first.second.HadMass();
    if ((m1+m2<MC) && (enforce || m1+m2<critM)) {
      wt = (sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * 
	    pow(sqr(m1+m2)/MC2,m_chi) * decit->second);
      if (wfweight>wfmax) {
	max1  = decit->first.first;
	max2  = decit->first.second;
	wfmax = wfweight;
      }
      totweight += wt;
    }
  }
  if (totweight<=0.) return 0.;
  had1 = had2 = Flavour(kf_none); 
  double disc(totweight * 0.9999999999*ran->Get());
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.HadMass();
    m2    = decit->first.second.HadMass();
    if ((m1+m2<MC) && (enforce || m1+m2<critM)) {
      wt = (sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * 
	    pow(sqr(m1+m2)/MC2,m_chi) * decit->second);
      disc -= wt;
      if (disc<0.) {
	had1 = decit->first.first;
	had2 = decit->first.second;
	break;
      }
    }
  }
  return wt/(16.*M_PI*MC*MC*MC);
}

void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Blob * blob,
				      const Flavour had1,const Flavour had2,
				      const bool & constrained)
{
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.HadMass()), m22 = sqr(had2.HadMass());

  cluster->BoostInCMSAndRotateOnZ();

  double E1((M2+m12-m22)/(2.*M)), pl2(sqr(E1)-m12);
  bool isbeam(false);
  double stheta, pt2;
  double masscor((cluster->GetTrip()->m_flav.HadMass() *
		  cluster->GetAnti()->m_flav.HadMass())/m_pt02);
  do { 
    stheta = 1.-2.*ran->Get(); 
    pt2    = pl2*sqr(stheta);
  } while (pt2>m_pt2max*m_pt2maxfac || 
	   sqr((*p_as)(pt2,false)/p_as->MaxValue())<ran->Get());
  double pt     = sqrt(pt2);
  int sign      = cluster->GetTrip()->m_mom[3]<0?-1:1;
  double pl1    = sign*sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi = cos(2.*M_PI*ran->Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1     = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2     = cluster->Momentum()-p1;

  if (p1[0]<0. || p2[0]<0.) throw Return_Value::Retry_Event;

  cluster->RotateAndBoostBack(p1);
  cluster->RotateAndBoostBack(p2);
  cluster->RotateAndBoostBack();

  Particle * left(new Particle(-1,had1,p1));
  left->SetNumber();
  left->SetInfo('P');
  left->SetFinalMass(had1.HadMass());
  Particle * right(new Particle(-1,had2,p2));
  right->SetNumber();
  right->SetInfo('P');
  right->SetFinalMass(had2.HadMass());
  control::s_AHAparticles+=2;

  if (blob!=NULL) {
    blob->AddToOutParticles(left);
    blob->AddToOutParticles(right);
  }
  // if (cluster->GetTrip()->m_info=='B' || cluster->GetAnti()->m_info=='B') {
  //   msg_Out()<<"==========================================================\n"
  // 	     <<"Cluster decay (pt = "<<pt<<") for cluster \n"<<(*cluster)<<"\n"
  // 	     <<"==> "<<left->Momentum()<<" + "<<right->Momentum()<<" for "
  // 	     <<left->Flav()<<" + "<<right->Flav()<<"\n"
  // 	     <<"==========================================================\n";
  // }
  if (m_ana) {
    Histogram* histo((m_histograms.find(std::string("PT_HH")))->second);
    histo->Insert(pt);
    Histogram* histo2((m_histograms.find(std::string("PT2_HH")))->second);
    histo2->Insert(pt*pt);
  }
}

bool Soft_Cluster_Handler::ShiftMomenta(Cluster_List * clin)
{
  if (!TryLocalCompensation(clin)) {
    return ForceMomenta(clin);
  }
  return true;
}

bool Soft_Cluster_Handler::TryLocalCompensation(Cluster_List * clin)
{
  int direx;
  Cluster * cluster, * partner;
  double mass1, mass2;
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    partner = NULL;
    direx   = 0;
    switch (cluster->size()) {
    case 1: 
      if (cluster->GetNBTrip()!=0 || cluster->GetNBAnti()!=0) {
	if (cluster->GetNBTrip()!=0 && cluster->GetNBAnti()!=0) {
	  if ((cluster->GetNBTrip()->size()!=1 && 
	       cluster->GetNBAnti()->size()!=1) ||
	      (cluster->GetNBTrip()->size()==1 && 
	       cluster->GetNBAnti()->size()==1)) {
	    if (1.>0.5) {
	      partner = cluster->GetNBAnti();
	      direx   = -1;
	    }
	    else {
	      partner = cluster->GetNBTrip();
	      direx   = +1;
	    }
	  }
	  else if (cluster->GetNBTrip()->size()==1 && 
		   cluster->GetNBAnti()->size()!=1) {
	    partner = cluster->GetNBTrip();
	    direx   = +1;
	  }
	  else {
	    partner = cluster->GetNBAnti();
	    direx   = -1;
	  }
	}
	else if (cluster->GetNBTrip()==0 && cluster->GetNBAnti()!=0) {
	  partner = cluster->GetNBAnti();
	}
	else if (cluster->GetNBTrip()!=0 && cluster->GetNBAnti()==0) {
	  partner = cluster->GetNBTrip();
	}
      }
      if (!partner) {
	return false;
      }
      mass1 = (*cluster)[0].HadMass();
      mass2 = partner->size()==1?(*partner)[0].HadMass():partner->Mass();
      if (sqr(mass1+mass2)>(cluster->Momentum()+partner->Momentum()).Abs2()) {
	if (direx==0) {
	  return false;
	}
	if (direx==-1) partner = cluster->GetNBTrip();
	if (direx== 1) partner = cluster->GetNBAnti();
	mass2 = partner->size()==1?(*partner)[0].HadMass():partner->Mass();
	if (sqr(mass1+mass2)>(cluster->Momentum()+partner->Momentum()).Abs2()) {
	  return false;
	}
      }
      masses.clear();
      momenta.clear();
      masses.push_back(mass1);
      masses.push_back(mass2);
      momenta.push_back(cluster->Momentum());
      momenta.push_back(partner->Momentum());
      if (!hadpars->AdjustMomenta(2,&momenta.front(),&masses.front())) {
	return false;
      }
      cluster->SetFlav((*cluster)[0]);
      cluster->SetMomentum(momenta[0]);
      partner->RescaleMomentum(momenta[1]);
      partner->SetMomentum(momenta[1]);
      break;
    case 2:
    case 0:
    default:
      break;
    }    
  }
  return true;
}

bool Soft_Cluster_Handler::ForceMomenta(Cluster_List * clin)
{
#ifdef AHAmomcheck
  Vec4D checkbef(SumMomentum(clin));
#endif
  size_t size(clin->size());
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  Cluster * cluster;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    switch (cluster->size()) {
    case 1: 
      masses.push_back((*cluster)[0].HadMass());
      break;
    case 2:
    case 0:
    default:
      masses.push_back(cluster->Mass());
      break;
    }
    momenta.push_back(cluster->Momentum());
  }

  if (!hadpars->AdjustMomenta(size,&momenta.front(),&masses.front())) {
    return false;
  }

  int pos(0);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) cluster->SetFlav((*cluster)[0]);
    else cluster->RescaleMomentum(momenta[pos]);
    cluster->SetMomentum(momenta[pos]);
    pos++;
  }
#ifdef AHAmomcheck
  Vec4D checkaft(SumMomentum(clin));
  double Q2(dabs((checkbef-checkaft).Abs2()));
  if (Q2>1.e-12 || IsNan(Q2)) {
    msg_Tracking()<<METHOD<<" yields a momentum violation for  "<<size<<" : \n"
		  <<"   "<<checkbef<<" - "<<checkaft<<" --> "
		  <<(checkbef-checkaft).Abs2()<<"("<<size<<").\n"
		  <<(*clin)<<"\n";
  }
  else msg_Debugging()<<METHOD<<" conserves momentum : "
		      <<(checkbef-checkaft).Abs2()<<"("<<size<<").\n";
#endif
  return true;
}

Vec4D Soft_Cluster_Handler::SumMomentum(Cluster_List * clin) {
  Cluster_Iterator cit;
  Vec4D listmom(0.,0.,0.,0.);
  for (cit=clin->begin();cit!=clin->end();cit++) listmom += (*cit)->Momentum(); 
  return listmom;
}

