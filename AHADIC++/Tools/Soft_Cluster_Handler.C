#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(bool ana) :
  p_as(hadpars->GetCoupling()),
  p_splitter(hadpars->GetSplitter()),
  p_singletransitions(hadpars->GetSingleTransitions()), 
  p_doubletransitions(hadpars->GetDoubleTransitions()),
  m_HHdecaymode(hadpars->Get(std::string("HHDecayMode"))),
  m_transitionoffset(hadpars->Get(std::string("Offset_C->H"))), 
  m_decayoffset(hadpars->Get(std::string("Offset_C->HH"))), 
  m_kappa(hadpars->Get(std::string("MassExponent_C->H"))), 
  m_lambda(hadpars->Get(std::string("WidthExponent_C->H"))), 
  m_chi(hadpars->Get(std::string("MassExponent_C->HH"))), 
  m_pt2max(sqr(hadpars->Get(string("ptmax")))),
  m_pt2maxfac(sqr(hadpars->Get(std::string("ptmax_factor")))),
  m_asform(p_as->Form()), m_pt02(p_as->PT02()), 
  m_transitions(0), m_dtransitions(0), m_decays(0), 
  m_forceddecays(0), m_lists(0),  
  m_ana(ana)
{
  switch (int(hadpars->Get(std::string("Selection_C->H")))) {
  case 3:
    m_tweightmode = Trans_Weight::waveonly;
    break;
  case 2:
    m_tweightmode = Trans_Weight::nonrelativistic;
    break;
  case 1:
  default:
    m_tweightmode = Trans_Weight::relativistic;
    break;
  }
  switch (int(hadpars->Get(std::string("Selection_C->HH")))) {
  case 14:
    m_dweightmode = Decay_Weight::waves;
    break;
  case 3:
    m_dweightmode = Decay_Weight::phasespace_masses_waves;
    break;
  case 2:
    m_dweightmode = Decay_Weight::phasespace_waves;
    break;
  case 1:
    m_dweightmode = Decay_Weight::phasespace;
    break;
  case 0:
  default:
    m_dweightmode = Decay_Weight::off;
  }

  if (m_ana) {
    m_histograms[string("PT_HH")]  = new Histogram(0,0.,10.,100);
    m_histograms[string("PT2_HH")] = new Histogram(0,0.,100.,2000);
    m_histograms[string("MassTransition")]       = new Histogram(0,0.,8.,100);
    m_histograms[string("HadronMassTransition")] = new Histogram(0,0.,8.,100);
  }
}

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{
  msg_Debugging()<<"@@@ "<<METHOD<<": "
		 <<m_transitions<<" transitions, "
		 <<m_dtransitions<<" double transitions, "
		 <<m_decays<<" decays and "
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
  if (clin->size()==1) return TreatSingleCluster(*clin->begin());
  if (!CheckListForTreatment(clin)) return true;
  m_lists++;
  double E(-1.);
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
    if (!EnforcedTransition(clin)) return false;
  }

  if (ShiftMomenta(clin)) {
    return AttachHadronsToBlob(clin,blob);
  }

  msg_Error()<<"Error in "<<METHOD<<" : \n"
	     <<"   Could not shift momenta.\n"
	     <<"   Will possibly lead to retrying the event.\n";
  return false;
}

bool Soft_Cluster_Handler::TreatSingleCluster(Cluster * cluster)
{
  switch (CheckCluster(cluster,true)) {
  case 1:
    msg_Debugging()<<"Potential problem in "<<METHOD<<" : \n"
		   <<"   Clusterlist with one element that "
		   <<"needs to transform to a hadron."
		   <<"\n"
		   <<"   Will possibly lead to retrying the event.\n";
    return false;
  case 2:
  case 0:
  default:
    break;
  }
  return true;
}

bool Soft_Cluster_Handler::TreatClusterDecay(Cluster_List * clin, Blob * blob)
{
  if (!CheckListForTreatment(clin)) return true;

  Cluster_Iterator cit=clin->begin();
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
    if (!EnforcedDecay(cluster,blob,true,clin)) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<clin->size()<<" clusters) : \n"
		 <<"   No enforced decay possible.\n";
      return false;
    }
    return true;
  }  
  if (UpdateClusterDecayKinematics((*clin->begin())->GetPrev())) {
    return AttachHadronsToBlob(clin,blob);
  }
  return false;
}

bool Soft_Cluster_Handler::
AttachHadronsToBlob(Cluster_List * clin,Blob * blob)
{
  Cluster_Iterator cit(clin->begin());
  Particle * part;
  Cluster * cluster;
  while (cit!=clin->end()) {
    cluster = (*cit);
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
    default:
      cit++;
      break;
    }
  }
  return true;
}

bool Soft_Cluster_Handler::UpdateClusterDecayKinematics(Cluster * cluster)
{
  Cluster * left(cluster->GetLeft()), * right(cluster->GetRight());
#ifdef AHAmomcheck
  cluster->CheckConsistency(msg_Error(),METHOD+std::string("entry"));
#endif
  cluster->BoostInCMSAndRotateOnZ();
  Vec4D  momleft(left->Momentum());
  double mass(cluster->Mass());
  double mass1((left->size()!=1)?left->Mass():(*left)[0].HadMass());
  double mass2((right->size()!=1)?right->Mass():(*right)[0].HadMass());
  if (left->size()==1 && right->size()==1) m_dtransitions += 1;

  double kt(sqrt(sqr(momleft[1])+sqr(momleft[2])));
  double cosphi(momleft[1]/kt), sinphi(momleft[2]/kt);
  double E1((sqr(mass)+sqr(mass1)-sqr(mass2))/(2.*mass)), E2(mass-E1);
  if (sqr(E1)-sqr(kt)-sqr(mass1)<0.) kt = sqrt(Max(0.,sqr(E1)-sqr(mass1)));
  if (sqr(E2)-sqr(kt)-sqr(mass2)<0.) kt = sqrt(Max(0.,sqr(E2)-sqr(mass2)));
  double pl1(sqrt(Max(0.,sqr(E1)-sqr(kt)-sqr(mass1))));
  double pl2(sqrt(Max(0.,sqr(E2)-sqr(kt)-sqr(mass2))));
  int sign(momleft[3]<0?-1:1);

  Vec4D momleft1 = Vec4D(E1,kt*cosphi,kt*sinphi,sign*pl1);
  Vec4D momright = Vec4D(E2,-kt*cosphi,-kt*sinphi,-sign*pl2); 
  if (left->size()==1)  left->SetMomentum(momleft1);
  else  left->RescaleMomentum(momleft1);
  if (right->size()==1) right->SetMomentum(momright);
  else right->RescaleMomentum(momright);

  cluster->RotateAndBoostBack();
#ifdef AHAmomcheck
  cluster->CheckConsistency(msg_Error(),METHOD+std::string("exit"));
#endif

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


int Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool lighter)
{
  msg_Tracking()<<METHOD<<"("<<cluster->Number()<<").\n";
  cluster->clear();

  Flavour haddec1(Flavour(kf_none)), haddec2(Flavour(kf_none));
  Flavour hadtrans(Flavour(kf_none));

  bool   direct(cluster->GetTrip()->m_info=='B' && 
		cluster->GetAnti()->m_info=='B');
  double decayweight(DecayWeight(cluster,haddec1,haddec2,direct));
  double transformweight(TransformWeight(cluster,hadtrans,lighter,false));
  if (decayweight>0. || transformweight>0.) {
    double totweight((Max(0.,decayweight)+Max(0.,transformweight))*0.9999999);
    if (transformweight<=0. || decayweight/totweight>ran->Get()) {
      m_decays      += 1;
      cluster->push_back(haddec1);
      cluster->push_back(haddec2);
      return 2;
    }
    if (transformweight>0.) {
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
  fl1 = fl1.Bar();
  fl2 = fl2.Bar();
  Proto_Particle *pp1(new Proto_Particle(fl1,cluster->GetTrip()->m_mom/2.,'l'));
  Proto_Particle *pp2(new Proto_Particle(fl2,cluster->GetTrip()->m_mom/2.,'l'));
  Proto_Particle *pp3(new Proto_Particle(fl3,cluster->GetAnti()->m_mom/2.,'l'));
  Proto_Particle *pp4(new Proto_Particle(fl4,cluster->GetAnti()->m_mom/2.,'l'));
  bool order(ran->Get()>0.5?true:false);
  Cluster cluster1((order?pp3:pp4),pp1), cluster2((!order?pp4:pp3),pp2);
  Flavour_Pair pair1, pair2;
  pair1.first = cluster1.GetTrip()->m_flav;
  pair1.first = cluster1.GetAnti()->m_flav;
  pair2.first = cluster2.GetTrip()->m_flav;
  pair2.first = cluster2.GetAnti()->m_flav;
  double mass(cluster->Mass());
  double wt1(TransformWeight(&cluster1,had1,false,true));
  double wt2(TransformWeight(&cluster2,had2,false,true));
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
	      Cluster_List * clin) {
  if (cluster->GetTrip()->m_info=='B' || cluster->GetAnti()->m_info=='B') {
    msg_Debugging()<<"@@@ "<<METHOD<<" for cluster "
		   <<"("<<cluster->GetTrip()->m_flav<<"/"
		   <<cluster->GetTrip()->m_info<<", "
		   <<cluster->GetAnti()->m_flav<<"/"
		   <<cluster->GetAnti()->m_info<<", "
		   <<"mass = "<<cluster->Mass()<<").\n";
  }
  Flavour had1,had2;
  double weight(DecayWeight(cluster,had1,had2,true)), weight1(-1.);
  if (weight<=0.) {
    if (cluster->GetTrip()->m_flav.IsDiQuark() && 
	cluster->GetAnti()->m_flav.IsDiQuark()) {
      if (!ClusterAnnihilation(cluster,had1,had2)) return false;
    }
    else {
      weight1 = TransformWeight(cluster,had1,true,true);
      if (weight1<=0.) {
	msg_Error()<<"Error in "<<METHOD<<" : \n"
		   <<"   No suitable single transition found, "
		   <<"will return false and hope for the best.\n";
	return false;
      }
      had2 = Flavour(kf_photon);
    }
    m_forceddecays++;m_decays--;
  }

  FixHHDecay(cluster,blob,had1,had2,constrained);
  m_decays++;
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
      msg_Error()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : \n"
		 <<"   Could not adjust momenta for : \n";
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Error()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Error()<<" ("<<((**cit)[0])<<" )";
	msg_Error()<<" for \n"<<(**cit)<<"\n";
	i++;
      }
      msg_Error()<<"   Will possibly lead to retrying the event.\n";
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
    msg_Error()<<METHOD<<" yields a momentum violation for  "<<size<<" : \n"
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
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;

  double MC(cluster->Mass());
  msg_Debugging()<<"@@ "<<METHOD<<" for cluster "<<cluster->Number()
		 <<" ("<<fpair.first<<", "
		 <<fpair.second<<", mass = "<<MC<<"): "<<endl
		 <<"@@ Heavy = "
		 <<p_singletransitions->GetHeaviestTransition(fpair)
		 <<" ("
		 <<(p_singletransitions->
		    GetHeaviestTransition(fpair).HadMass())<<"), "
		 <<"light = "
		 <<p_singletransitions->GetLightestTransition(fpair)
		 <<" ("
		 <<(p_singletransitions->
		    GetLightestTransition(fpair).HadMass())<<"), "
		 <<"enforce = "<<enforce<<", lighter = "
		 <<lighter<<" for "<<hadron<<".\n";

  if (!enforce && 
      p_singletransitions->GetHeaviestMass(fpair)<MC+m_transitionoffset) {
    hadron = Flavour(kf_none);
    return 0.;
  }
  if ((enforce || lighter) && p_singletransitions->GetLightestMass(fpair)>MC) {
    msg_Error()<<"Error in "<<METHOD<<"("<<lighter<<", "<<enforce<<") :\n"
	       <<"   Cluster too light, no transformation possible.\n"
	       <<"   Cluster = "<<cluster->Number()<<" ["
	       <<fpair.first<<"("<<fpair.first.HadMass()<<"), "
	       <<fpair.second<<"("<<fpair.second.HadMass()<<"), "
	       <<"mass = "<<MC<<"] vs. "
	       <<"lightest mass = "
	       <<p_singletransitions->GetLightestMass(fpair)<<" for "
	       <<p_singletransitions->GetLightestTransition(fpair)<<".\n"
	       <<(*cluster)<<"\n";
    return 0.;
  }
  Single_Transition_Miter stiter = 
    p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential error in  "<<METHOD<<" :"<<endl
	       <<"   Did not find any entry for "
	       <<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<"\n"
	       <<"   will continue and hope for the best."<<endl;
    hadron = Flavour(kf_none);
    return -1.;
  }
  Single_Transition_List * stl(stiter->second);
  Single_Transition_Siter  start(stl->begin()),siter;
  if (lighter) {
    if (hadron!=Flavour(kf_none)) {
      do {
	if (start->first==hadron) {
	  siter = start;
	  siter--;
	  if ((++siter)!=stl->end()) start++;
	  else return 0.;
	  break;
	}
	else start++;
      } while (start!=stl->end());
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

  if (lighter || 
      cluster->GetTrip()->m_info=='B' || 
      cluster->GetAnti()->m_info=='B') {
    msg_Debugging()<<"@@     --> "<<hadron<<" ("<<hadron.HadMass()<<", "
		   <<"wt = "<<disc/totweight<<").\n";
  }
  return wt/(16.*M_PI*MC);
}

double Soft_Cluster_Handler::
TransformKin(const double MC,const Flavour & flav,
	     const bool & enforce) {
  double weight(1.);
  double mass(flav.HadMass()),mass2(mass*mass);
  double width(Max(flav.Width(),1.e-6)), width2(width*width);
  switch (m_tweightmode) {
  case Trans_Weight::waveonly:
    break;
  case Trans_Weight::nonrelativistic:
    if (!enforce && dabs(MC-mass)>20.*width) return 0.;
    weight *= pow(mass2/(sqr(MC-mass) + width2/4.),m_kappa);
    weight *= pow(mass*width/(sqr(MC-mass) + width2/4.),m_lambda);
    break;
  case Trans_Weight::relativistic:
  default:
    if (!enforce && dabs(MC-mass)>20.*width) return 0.;
    weight *= pow(sqr(mass2)/(sqr(MC*MC-mass2) + mass2*width2),m_kappa);
    weight *= pow(mass2*width2/(sqr(MC*MC-mass2) + mass2*width2),m_lambda);
    break;
  }
  return weight;
}

double Soft_Cluster_Handler::
DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2,
	    const bool & enforce)
{
  if (m_dweightmode==Decay_Weight::off && !enforce) return 0.;
  Flavour_Pair flpair;
  flpair.first  = cluster->GetTrip()->m_flav;
  flpair.second = cluster->GetAnti()->m_flav;

  Double_Transition_Miter dtliter = 
    p_doubletransitions->GetTransitions()->find(flpair);
  double MC(cluster->Mass());
  
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential Error in "<<METHOD<<" : "<<endl
	       <<"   No viable transition found for "
	       <<flpair.first<<"/"<<flpair.second
	       <<" (mass = "<<MC<<") in "
	       <<p_doubletransitions->GetTransitions()->size()<<"."<<endl
	       <<"   Return 'false' and hope for the best.\n";
    return 0.;
  }
  if (p_doubletransitions->GetLightestMass(flpair)>MC) {
    if (enforce) {
      msg_Tracking()<<"Warning in "<<METHOD<<" : "<<endl
		    <<"   No viable transition found for "
		    <<flpair.first<<"/"<<flpair.second
		    <<" (mass = "<<MC<<")."<<endl
		    <<"   Return 'false' and hope for the best.\n";
    }
    return 0.;
  }
  if (!enforce && 
      p_doubletransitions->GetLightestMass(flpair)*(1.-m_decayoffset)+
      p_doubletransitions->GetHeaviestMass(flpair)*m_decayoffset<MC) {
    msg_Debugging()<<"@@      --> too heavy, no decay.\n";
    had1 = had2 = Flavour(kf_none);
    return 0.;
  }
  else {
    msg_Tracking()<<"@@ "<<METHOD<<" for cluster "<<cluster->Number()
		  <<" ("<<flpair.first<<", "<<flpair.second<<", "
		  <<"mass = "<<MC<<"): "<<endl
		  <<"@@ Heavy = "
		  <<p_doubletransitions->GetHeaviestTransition(flpair).first
		  <<" + "
		  <<p_doubletransitions->GetHeaviestTransition(flpair).second
		  <<" ("<<(p_doubletransitions->
			   GetHeaviestTransition(flpair).first.HadMass()+
			   p_doubletransitions->
			   GetHeaviestTransition(flpair).second.HadMass())
		  <<"), "
		  <<"light = "
		  <<p_doubletransitions->GetLightestTransition(flpair).first
		  <<" + "
		  <<p_doubletransitions->GetLightestTransition(flpair).second
		  <<" ("
		  <<(p_doubletransitions->
		     GetLightestTransition(flpair).first.HadMass()+
		     p_doubletransitions->
		     GetLightestTransition(flpair).second.HadMass())<<"), "
		  <<"enforce = "<<enforce<<".\n";
  }

  double  totweight(0.),MC2(MC*MC),m1,m2,wt(1.),wfweight(0.),wfmax(0.);
  Flavour max1, max2;
  double  expo1(m_chi),expo2(m_chi);
  //if (!enforce &&
  //    (cluster->GetTrip()->m_info=='B' || cluster->GetTrip()->m_info=='L'))
  //  expo1 = 0.;
  //if (!enforce &&
  //    (cluster->GetAnti()->m_info=='B' || cluster->GetAnti()->m_info=='L'))
  //  expo2 = 0.;
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.HadMass();
    m2  = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = 1.;
      if (m_dweightmode==Decay_Weight::phasespace ||
	  m_dweightmode==Decay_Weight::phasespace_waves ||
	  m_dweightmode==Decay_Weight::phasespace_masses_waves)
	wt *= sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
      if (m_dweightmode==Decay_Weight::phasespace_masses_waves)
	wt *= pow(2.*m1/MC,expo1)*pow(2.*m2/MC,expo2);
      if (m_dweightmode==Decay_Weight::phasespace_masses_waves ||
	  m_dweightmode==Decay_Weight::phasespace_waves ||
	  m_dweightmode==Decay_Weight::waves)
	wt *= wfweight = decit->second;
      if (wfweight>wfmax) {
	max1  = decit->first.first;
	max2  = decit->first.second;
	wfmax = wfweight;
      }
      totweight += wt;
    }
  }
  if (totweight<=0.) {
    msg_Debugging()<<"Error in "<<METHOD<<" :\n"
		   <<"   Cluster of mass "<<MC<<" from {"<<flpair.first<<", "
		   <<flpair.second<<"} passed mass conditions,\n"
		   <<"   but no viable transition found.\n"
		   <<"   Return 0 and hope for the best.\n";
    return 0.;
  }

  had1 = had2 = Flavour(kf_none); 
  double disc(totweight * 0.9999999999*ran->Get());
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.HadMass();
    m2    = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = 1.;
      if (m_dweightmode==Decay_Weight::phasespace ||
	  m_dweightmode==Decay_Weight::phasespace_waves ||
	  m_dweightmode==Decay_Weight::phasespace_masses_waves)
	wt *= sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
      if (m_dweightmode==Decay_Weight::phasespace_masses_waves)
	wt *= pow(2.*m1/MC,expo1)*pow(2.*m2/MC,expo2);
      if (m_dweightmode==Decay_Weight::phasespace_masses_waves ||
	  m_dweightmode==Decay_Weight::phasespace_waves ||
	  m_dweightmode==Decay_Weight::waves)
	wt *= wfweight = decit->second;
      disc -= wt;
      if (disc<0.) {
	had1 = decit->first.first;
	had2 = decit->first.second;
	break;
      }
    }
  }
  msg_Tracking()<<"@@ ["<<cluster->Number()<<"]  --> "<<had1<<" + "<<had2<<" "
		<<"( waves = "<<wfweight<<", max = "<<wfmax
		<<" from "<<max1<<" + "<<max2<<"), "
		<<"expos = {"<<expo1<<", "<<expo2<<"}.\n";
  return wt/(16.*M_PI*MC*MC*MC);
}

void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Blob * blob,
				      const Flavour had1,const Flavour had2,
				      const bool & constrained)
{
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.HadMass()), m22 = sqr(had2.HadMass());

  cluster->BoostInCMSAndRotateOnZ();

  double E1((M2+m12-m22)/(2.*M)), p2max(sqr(E1)-m12);
  double pt2max(p2max);
  bool isbeam(false);
  if (cluster->GetTrip()->m_info=='B' || 
      cluster->GetAnti()->m_info=='B') {  
    pt2max = Min(pt2max,Min(cluster->GetTrip()->m_kt2max,
			    cluster->GetAnti()->m_kt2max));
    isbeam = true;
  }
  else if (cluster->GetTrip()->m_info=='L' || 
	   cluster->GetAnti()->m_info=='L') {
    pt2max = Min(pt2max,
		 Min(cluster->GetTrip()->m_kt2max*
		     m_pt02/Max(m_pt02,sqr(cluster->GetTrip()->m_mass)),
		     cluster->GetAnti()->m_kt2max*
		     m_pt02/Max(m_pt02,sqr(cluster->GetAnti()->m_mass))));
  }
  else if (cluster->GetTrip()->m_info!='L' || 
	   cluster->GetAnti()->m_info!='L') {
    pt2max = Min(pt2max,
		 Min(cluster->GetTrip()->m_kt2max*
		     m_pt02/Max(m_pt02,sqr(cluster->GetTrip()->m_mass)/4.),
		     cluster->GetAnti()->m_kt2max*
		     m_pt02/Max(m_pt02,sqr(cluster->GetAnti()->m_mass)/4.)));
  }

  
  if (IsNan(pt2max) || pt2max<0.) 
    msg_Error()<<"Error in "<<METHOD
	       <<"(pt2max = "<<pt2max<<", p2max = "<<p2max<<") for "
	       <<M<<" --> "<<sqrt(m12)<<"("<<had1<<") + "
	       <<sqrt(m22)<<"("<<had2<<"),\n"<<(*cluster)<<".\n";

  double pt(0.),pl1(0.);
  if (m_HHdecaymode==0) {
    double ctheta = 1.-2.*ran->Get(), stheta=sqrt(1.-ctheta*ctheta);
    pt            = sqrt(pt2max)*stheta;
    pl1           = p2max*ctheta;
  }
  else {
    double masscor = Max(sqr(cluster->GetTrip()->m_mass)/m_pt2max,1.) *
      Max(sqr(cluster->GetAnti()->m_mass)/m_pt2max,1.);
    double pt2     = SelectPT2(pt2max,masscor,true,isbeam);
    pt             = sqrt(pt2);
    int sign       = cluster->GetTrip()->m_mom[3]<0?-1:1;
    pl1            = sign*sqrt(sqr(E1)-sqr(pt)-m12);
  }
  double cosphi  = cos(2.*M_PI*ran->Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  if (p1[0]<0. || p2[0]<0.) {
    msg_Error()<<"Error in "<<METHOD<<": negative hadron energies\n"
	       <<(*cluster)<<"\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }

  if (cluster->GetLeft()) {
    delete cluster->GetLeft()->GetTrip();
    delete cluster->GetLeft()->GetAnti();
    delete cluster->GetLeft();
  }
  if (cluster->GetRight()) {
    delete cluster->GetRight()->GetTrip();
    delete cluster->GetRight()->GetAnti();
    delete cluster->GetRight();
  }

  cluster->SetLeft(new Cluster(p1,had1,false));
  cluster->GetLeft()->SetPrev(cluster);
  cluster->SetRight(new Cluster(p2,had2,false));
  cluster->GetRight()->SetPrev(cluster);
  
  cluster->RotateAndBoostBack();

  Particle * left(cluster->GetLeft()->GetSelf());
  left->SetFinalMass();
  Particle * right(cluster->GetRight()->GetSelf());
  right->SetFinalMass();
  //  if (show) {
  msg_Tracking()<<METHOD<<"("<<m_HHdecaymode<<"): "
		<<"attach two hadrons for cluster ("
		<<cluster->Number()<<"; {"<<cluster->GetTrip()->m_info
		<<cluster->GetAnti()->m_info<<"}), "
		<<pt<<" from {"<<sqrt(pt2max)<<", "<<sqrt(m_pt2max)<<"}, "
		<<cluster->Momentum()
		<<",mass = "<<cluster->Mass()<<".\n"
		<<"    ("<<left->Flav()<<", "<<left->Momentum()
		<<" --> pt = "<<left->Momentum().PPerp()<<", "
		<<"y = "<<left->Momentum().Y()<<")\n"
		<<"    ("<<right->Flav()<<", "<<right->Momentum()
		<<" --> pt = "<<right->Momentum().PPerp()<<", "
		<<"y = "<<right->Momentum().Y()<<").\n";
  //  }

  if (blob!=NULL) {
    blob->AddToOutParticles(left);
    blob->AddToOutParticles(right);
    delete cluster->GetLeft();
    delete cluster->GetRight();
  }

  if (m_ana) {
    Histogram* histo((m_histograms.find(std::string("PT_HH")))->second);
    histo->Insert(pt);
    Histogram* histo2((m_histograms.find(std::string("PT2_HH")))->second);
    histo2->Insert(pt*pt);
  }
}

double Soft_Cluster_Handler::
SelectPT2(const double & pt2max,const double & masscor,const bool & expo,const bool & beam) {
  double pt2(0.), arg((pt2max+m_pt02)/m_pt02), exparg(m_pt2maxfac*m_pt2max/masscor), wt(0.);
  bool   runit(true);
  do {
    pt2 = m_pt02*(pow(arg,ran->Get())-1.);
    wt  = p_as->Weight(pt2);
    if (expo && pt2>m_pt2max) 
      wt *= exp(-sqr((pt2-m_pt2max)/exparg));
  } while (wt<ran->Get());
  return pt2;
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
    if (size>1) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : \n"
		 <<"   Could not adjust momenta for : \n";
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Error()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Error()<<" ("<<((**cit)[0])<<" )";
	msg_Error()<<" for \n"<<(**cit)<<"\n";
	i++;
      }
      msg_Error()<<"   Will possibly lead to retrying the event.\n";
    }
    return false;
  }

  int pos(0);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      cluster->SetFlav((*cluster)[0]);
    }
    else {
      cluster->RescaleMomentum(momenta[pos]);
    }
    cluster->SetMomentum(momenta[pos]);
    pos++;
  }
#ifdef AHAmomcheck
  Vec4D checkaft(SumMomentum(clin));
  double Q2(dabs((checkbef-checkaft).Abs2()));
  if (Q2>1.e-12 || IsNan(Q2)) {
    msg_Error()<<METHOD<<" yields a momentum violation for  "<<size<<" : \n"
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

