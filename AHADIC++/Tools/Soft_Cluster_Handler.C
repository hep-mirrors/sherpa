#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(bool ana) :
  p_as(hadpars.GetCoupling()),
  p_splitter(hadpars.GetSplitter()),
  p_singletransitions(hadpars.GetSingleTransitions()), 
  p_doubletransitions(hadpars.GetDoubleTransitions()),
  m_transitionoffset(hadpars.Get(std::string("Offset_C->H"))), 
  m_decayoffset(hadpars.Get(std::string("Offset_C->HH"))), 
  m_kappa(hadpars.Get(std::string("MassExponent_C->H"))), 
  m_lambda(hadpars.Get(std::string("WidthExponent_C->H"))), 
  m_chi(hadpars.Get(std::string("MassExponent_C->HH"))), 
  m_pt2min(dabs(hadpars.Get(string("pt2min")))),
  m_pt2maxfac(sqr(hadpars.Get(std::string("ptmax_factor")))),
  m_transitions(0), m_dtransitions(0), m_decays(0), m_forceddecays(0), m_lists(0),  
  m_ana(ana)
{
  switch (int(hadpars.Get(std::string("Selection_C->H")))) {
  case 3:
    m_tweightmode = TransWeight::waveonly;
    break;
  case 2:
    m_tweightmode = TransWeight::nonrelativistic;
    break;
  case 1:
  default:
    m_tweightmode = TransWeight::relativistic;
    break;
  }
  switch (int(hadpars.Get(std::string("Selection_C->HH")))) {
  case 14:
    m_dweightmode = DecayWeight::waves;
    break;
  case 3:
    m_dweightmode = DecayWeight::phasespace_masses_waves;
    break;
  case 2:
    m_dweightmode = DecayWeight::phasespace_waves;
    break;
  case 1:
    m_dweightmode = DecayWeight::phasespace;
    break;
  case 0:
  default:
    m_dweightmode = DecayWeight::off;
  }

  if (m_ana) {
    m_histograms[string("PT_HH")]                = new Histogram(0,0.,2.0,100);
    m_histograms[string("PT2_HH")]               = new Histogram(0,0.,4.0,200);
    m_histograms[string("MassTransition")]       = new Histogram(0,0.,8.,100);
    m_histograms[string("HadronMassTransition")] = new Histogram(0,0.,8.,100);
  }
}

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{
  msg_Out()<<"@@@ "<<METHOD<<": "
	   <<m_transitions<<" transitions, "<<m_dtransitions<<" double transitions, "
	   <<m_decays<<" decays and "<<m_forceddecays<<" forced decays."<<std::endl
	   <<"@@@ "<<METHOD<<": "<<m_lists<<" transitions from original dipole list."<<std::endl;
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

bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin, Blob * blob,bool cludec)
{
  msg_Tracking()<<std::endl<<"@@@@ "<<METHOD<<" with "<<clin->size()<<" clusters:"<<std::endl
		<<(*clin)<<std::endl;
  if (clin->size()==1)              return TreatSingleCluster(*clin->begin(),blob);
  if (cludec && clin->size()==2)    return TreatClusterDecay(clin,blob);

  if (!CheckListForTreatment(clin)) return true;
  m_lists++;
  double E(-1.);

  msg_Tracking()<<"@@@@ "<<METHOD<<" for "<<clin->size()<<" clusters to deal with."<<std::endl;
  while (!CheckIfAllowed(clin,E)) {
    if (UpdateTransitions(clin)) continue;
    if (EnforcedTransition(clin)) break;
    msg_Debugging()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Cannot treat list, this may trigger a new event."<<std::endl;
    return false;
  }
  if (ShiftMomenta(clin)) {
    //msg_Out()<<(*clin)<<std::endl;
    AttachHadronsToBlob(clin,blob);
    //msg_Out()<<"@@@ "<<METHOD<<": Fragmentation blob now reads: "<<std::endl<<(*blob)<<std::endl;
    return true;
  }

  msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	     <<"   Could not shift momenta."<<std::endl
	     <<"   Will possibly lead to retrying the event."<<std::endl;
  return false;
}

bool Soft_Cluster_Handler::TreatSingleCluster(Cluster * cluster, Blob * blob)
{
  msg_Tracking()<<std::endl<<"############# "<<METHOD<<" for "<<cluster->Number()<<"."<<std::endl;
  switch (CheckCluster(cluster,true)) {
  case -1:
  case 1:
    msg_Tracking()<<"Potential problem in "<<METHOD<<" : "<<std::endl
		  <<"   Clusterlist with one element that needs to transform to a hadron."
		  <<std::endl
		  <<"   Will possibly lead to retrying the event."<<std::endl;
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
  msg_Tracking()<<std::endl<<"@@@@@@ "<<METHOD<<" for "<<clin->size()<<" clusters."<<std::endl;
  if (!CheckListForTreatment(clin)) return true;

  double E(-1.);
  while (!CheckIfAllowed(clin,E)) {
    if (UpdateTransitions(clin)) continue;
    Cluster_Iterator cit=clin->begin();
    Cluster * left((*cit++)), * right((*cit));
    msg_Tracking()<<"@@@@@@ "<<METHOD<<": check this "<<left<<"/"<<right<<"."<<std::endl;
    if (left->GetPrev()!=right->GetPrev()) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<clin->size()<<" clusters) : "<<std::endl
		 <<"   No common previous cluster."<<std::endl;
      return false;
    }
    Cluster * cluster(left->GetPrev());
    msg_Tracking()<<"@@@@@ clear cluster list: "<<clin->size()<<" elements."<<std::endl;
    clin->clear();
    msg_Tracking()<<"@@@@@ now in cluster list: "<<clin->size()<<" elements."<<std::endl;
    if (!EnforcedDecay(cluster,blob,clin)) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<clin->size()<<" clusters) : "<<std::endl
		 <<"   No enforced decay possible."<<std::endl;
      return false;
    }
    msg_Tracking()<<"@@@@@@ "<<METHOD<<" enforced decay worked out."<<std::endl;
    return true;
  }  
  if (UpdateClusterDecayKinematics((*clin->begin())->GetPrev()))
    return AttachHadronsToBlob(clin,blob);
  return false;
}

bool Soft_Cluster_Handler::AttachHadronsToBlob(Cluster_List * clin,ATOOLS::Blob * blob)
{
  msg_Tracking()<<"$$$$$$ "<<METHOD<<" for "<<clin->size()<<" clusters."<<std::endl;
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
      msg_Tracking()<<"$$ attach one hadron ("<<part->Flav()<<") from cluster "<<cluster->Number()<<"."<<std::endl;
      delete cluster->GetTrip();
      delete cluster->GetAnti();
      delete cluster;
      cit = clin->erase(cit);
      break;
    case 2:
      msg_Tracking()<<"$$ will attach two hadrons ("<<(*cluster)[0]<<" + "<<(*cluster)[1]<<") "
		    <<"from cluster "<<cluster->Number()<<"."<<std::endl;
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
  msg_Tracking()<<"@@ in "<<METHOD<<":"<<std::endl;
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
  msg_Tracking()<<"@@  Original momenta = "
		<<left->Momentum()<<" & "<<right->Momentum()<<"."<<std::endl
		<<"@@  Masses "<<mass<<" -> "
		<<mass1<<"("<<(left->size()==1)<<") + "
		<<mass2<<"("<<(right->size()==1)<<")."<<std::endl;
  double kt(sqrt(sqr(momleft[1])+sqr(momleft[2])));
  double cosphi(momleft[1]/kt), sinphi(momleft[2]/kt);
  double E1((sqr(mass)+sqr(mass1)-sqr(mass2))/(2.*mass)), E2(mass-E1);
  if (sqr(E1)-sqr(kt)-sqr(mass1)<0.) kt = sqrt(Max(0.,sqr(E1)-sqr(mass1)));
  if (sqr(E2)-sqr(kt)-sqr(mass2)<0.) kt = sqrt(Max(0.,sqr(E2)-sqr(mass2)));
  double pl1(sqrt(Max(0.,sqr(E1)-sqr(kt)-sqr(mass1))));
  double pl2(sqrt(Max(0.,sqr(E2)-sqr(kt)-sqr(mass2))));
  msg_Tracking()<<"@@  Energies = "<<E1<<" & "<<E2<<", plongs = "<<pl1<<" & "<<pl2<<"."<<std::endl;
  int sign(momleft[3]<0?-1:1);

  Vec4D momleft1 = Vec4D(E1,kt*cosphi,kt*sinphi,sign*pl1);
  Vec4D momright = Vec4D(E2,-kt*cosphi,-kt*sinphi,-sign*pl2); 
  msg_Tracking()<<"@@  New momenta      = "<<left->Momentum()<<" & "<<right->Momentum()<<"."<<std::endl;
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
  msg_Tracking()<<"############# "<<METHOD<<" for "<<clin->size()<<" clusters."<<std::endl;
  Cluster_Iterator cit;
  Cluster * cluster;
  int    size(0),count(0);
  for (cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster==NULL || !cluster->Active()) continue;
    size = CheckCluster(cluster,false);
    switch (size) {
    case -1:
      return false;
    default:
      count += size;
    }
  }
  if (count==0) return false;
  return true;
}


int Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool lighter)
{
  cluster->clear();

  Flavour haddec1(Flavour(kf_none)), haddec2(Flavour(kf_none)), hadtrans(Flavour(kf_none));
  msg_Tracking()<<"@@@@ "<<METHOD<<" for cluster ("<<cluster->Number()<<":"
  		<<cluster->GetTrip()->m_flav<<" "<<cluster->GetAnti()->m_flav<<", "
  		<<"m = "<<cluster->Mass()<<")"<<std::endl;
  double decayweight(DecayWeight(cluster,haddec1,haddec2,false));
  if (decayweight>0.) 
    msg_Tracking()<<"@@@@  --> "<<haddec1<<" + "<<haddec2<<" with weight = "<<decayweight<<"."<<std::endl;
  double transformweight(TransformWeight(cluster,hadtrans,lighter,false));
  if (transformweight>0.) 
    msg_Tracking()<<"@@@@  --> "<<hadtrans<<" with weight = "<<transformweight<<"."<<std::endl;
  if (decayweight>0. || transformweight>0.) {
    double totweight((decayweight+transformweight)*0.9999999);
    if (totweight<=0. || decayweight/totweight>ran.Get()) {
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
  msg_Tracking()<<"@@@@  --> NONE with weight = "<<transformweight<<"."<<std::endl;  
  return 0;
}

bool Soft_Cluster_Handler::CheckIfAllowed(Cluster_List * clin,double & E) {
  double totmass(0.);
  Vec4D  totmom(0.,0.,0.,0.);
  Cluster * cluster;
  std::vector<Flavour> flavs;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    switch (cluster->size()) {
    case 1: 
      totmass += (*cluster)[0].HadMass();
      flavs.push_back((*cluster)[0]);
      break;
    case 2:
    case 0:
    default:
      flavs.push_back(Flavour(kf_cluster));
      totmass += cluster->Mass();
      break;
    }
    if (E<0) totmom += cluster->Momentum();
  }
  if (E<0) E = sqrt(totmom.Abs2());
  if (clin->size()==2) {
    if (totmass<E) 
      msg_Tracking()<<"++ "<<METHOD<<" : transitions to "<<flavs[0]<<" & "<<flavs[1]<<" allowed."<<std::endl;
    else
      msg_Tracking()<<"++ "<<METHOD<<" : transitions to "<<flavs[0]<<" & "<<flavs[1]<<" not allowed."<<std::endl;
  }
  return (totmass<E);
}

bool Soft_Cluster_Handler::UpdateTransitions(Cluster_List * clin) {
  msg_Tracking()<<"@@@@ "<<METHOD<<" for "<<clin->size()<<" clusters. "<<std::endl;
  Cluster * cluster, * winner(NULL);
  Flavour hadron,winhad;
  double  tfwt, maxwt(0.);
  bool    found(false);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      hadron = (*cluster)[0];
      tfwt   = TransformWeight(cluster,hadron,true,false);
      if (tfwt>maxwt) {
	winner = cluster;
	winhad = hadron;
	maxwt  = tfwt;
	found  = true; 
      }
    }
  }
  if (found) {
    msg_Tracking()<<"@@@@ "<<METHOD<<" found a winner: "<<winner<<"."<<std::endl;
    msg_Tracking()<<"@@@@ Replace "<<(*winner)[0]<<" with "<<winhad<<"."<<std::endl;
    winner->clear();
    winner->push_back(winhad);
    return true;
  }
  msg_Tracking()<<"@@@@ "<<METHOD<<": no winner found."<<std::endl;
  return false;
}

bool Soft_Cluster_Handler::ClusterAnnihilation(Cluster * cluster,Flavour & had1,Flavour & had2) {
  int kfc1(int(cluster->GetTrip()->m_flav.Kfcode())); 
  int kfc2(int(cluster->GetAnti()->m_flav.Kfcode())); 
  kf_code kfc11(kfc1/1000),kfc12((kfc1-kfc11*1000)/100);
  kf_code kfc21(kfc2/1000),kfc22((kfc2-kfc21*1000)/100);
  Flavour fl1(kfc11), fl2(kfc12), fl3(kfc21), fl4(kfc22);
  fl1 = fl1.Bar();
  fl2 = fl2.Bar();
  Proto_Particle * pp1 = new Proto_Particle(fl1,cluster->GetTrip()->m_mom/2.,'l');
  Proto_Particle * pp2 = new Proto_Particle(fl2,cluster->GetTrip()->m_mom/2.,'l');
  Proto_Particle * pp3 = new Proto_Particle(fl3,cluster->GetAnti()->m_mom/2.,'l');
  Proto_Particle * pp4 = new Proto_Particle(fl4,cluster->GetAnti()->m_mom/2.,'l');
  bool order(ran.Get()>0.5?true:false);
  Cluster cluster1((order?pp3:pp4),pp1), cluster2((!order?pp4:pp3),pp2);
  Flavour_Pair pair1, pair2;
  pair1.first = cluster1.GetTrip()->m_flav;
  pair1.first = cluster1.GetAnti()->m_flav;
  pair2.first = cluster2.GetTrip()->m_flav;
  pair2.first = cluster2.GetAnti()->m_flav;
  double mass(cluster->Mass());
  double wt1(TransformWeight(&cluster1,had1,false,true));
  double wt2(TransformWeight(&cluster2,had2,false,true));
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
  //std::cout<<"====================================================="<<std::endl
  //	   <<" "<<METHOD<<" for "
  //	   <<cluster->GetTrip()->m_flav<<"/"<<cluster->GetAnti()->m_flav
  //	   <<" with mass = "<<mass<<" --> "<<std::endl
  //	   <<" "<<kfc1<<"/"<<kfc2<<" --> "
  //	   <<kfc11<<"+"<<kfc12<<"/"<<kfc21<<"+"<<kfc22<<"  --> "
  //	   <<fl1<<"+"<<fl2<<"/"<<fl3<<"+"<<fl4<<std::endl
  //	   <<" --> clusters : "<<std::endl<<cluster1<<cluster2
  //	   <<" --> hadrons : "<<had1<<" / "<<had2<<std::endl
  //	   <<"====================================================="<<std::endl;
  return true;
}

bool Soft_Cluster_Handler::EnforcedDecay(Cluster * cluster, Blob * blob,
					 Cluster_List * clin) {
  msg_Tracking()<<"@@@ "<<METHOD<<" for cluster "
		<<"("<<cluster->GetTrip()->m_flav<<"/"<<cluster->GetAnti()->m_flav<<", "
		<<"mass = "<<cluster->Mass()<<")."<<std::endl;
  Flavour had1,had2;
  double weight(DecayWeight(cluster,had1,had2,true)), weight1(-1.);
  if (weight<=0.) {
    if (cluster->GetTrip()->m_flav.IsDiQuark() && cluster->GetAnti()->m_flav.IsDiQuark()) {
      if (!ClusterAnnihilation(cluster,had1,had2)) {
	// 	std::cout<<"====================================================="<<std::endl
	// 		 <<" "<<METHOD<<" for "
	// 		 <<cluster->GetTrip()->m_flav<<"/"<<cluster->GetAnti()->m_flav
	// 		 <<" with mass = "<<cluster->Mass()<<"."<<std::endl
	// 		 <<"====================================================="<<std::endl;
	return false;
      }
      //       std::cout<<"====================================================="<<std::endl
      // 	       <<" "<<METHOD<<" continues with "
      // 	       <<cluster->GetTrip()->m_flav<<"/"<<cluster->GetAnti()->m_flav
      // 	       <<", mass = "<<cluster->Mass()<<" --> "<<had1<<" & "<<had2<<"."<<std::endl
      // 	       <<"====================================================="<<std::endl;
    }
    else {
      weight1 = TransformWeight(cluster,had1,true,true);
      if (weight1<=0.) {
	msg_Tracking()<<"@@@ "<<METHOD<<" no viable decay, no viable transformation."<<std::endl;
	return false;
      }
      had2 = Flavour(kf_photon);
    }
    m_forceddecays++;m_decays--;
  }
  msg_Tracking()<<"@@@ "<<METHOD<<" enforced decay to "<<had1<<"/"<<had2<<" (weights = "
		<<weight<<", "<<weight1<<")."<<std::endl;
  FixHHDecay(cluster,blob,had1,had2);
  m_decays++;
  return true;
}

bool Soft_Cluster_Handler::EnforcedTransition(Cluster_List * clin) {
  bool success(true);
#ifdef AHAmomcheck
  Vec4D checkbef(SumMomentum(clin));
#endif
  size_t size(0);
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  Cluster * cluster;
  //std::cout<<"{{{ "<<METHOD<<" : }}}"<<std::endl;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    //std::cout<<"   "<<cluster->Number()<<" ("<<cluster->size()<<") ";
    switch (cluster->size()) {
    case 1: 
      //std::cout<<"with "<<(*cluster)[0]<<"."<<std::endl;
      masses.push_back((*cluster)[0].HadMass());
      momenta.push_back(cluster->Momentum());
      ++size;
      break;
    case 2:
    case 0:
      //std::cout<<"."<<std::endl;
      masses.push_back(sqrt(Max(cluster->GetTrip()->m_mom.Abs2(),0.)));
      masses.push_back(sqrt(Max(cluster->GetAnti()->m_mom.Abs2(),0.)));
      momenta.push_back(cluster->GetTrip()->m_mom);
      momenta.push_back(cluster->GetAnti()->m_mom);
      size+=2;
    default:
      break;
    }
  }
  if (!hadpars.AdjustMomenta(size,&momenta.front(),&masses.front())) {
    if (size>1 /*&& msg->LevelIsDebugging()*/) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : "<<std::endl
		    <<"   Could not adjust momenta for : "<<std::endl;
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Error()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Error()<<" ("<<((**cit)[0])<<" )";
	msg_Error()<<" for "<<std::endl<<(**cit)<<std::endl;
	i++;
      }
      msg_Error()<<"   Will possibly lead to retrying the event."<<std::endl;
    }
    exit(1);
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
      if (cluster->Mass()<(*cluster)[0].HadMass()+(*cluster)[1].HadMass()) success=false;
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
    msg_Error()<<METHOD<<" yields a momentum violation for  "<<size<<" : "<<std::endl
	       <<"   "<<checkbef<<" - "<<checkaft<<" --> "
	       <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl
	       <<(*clin)<<std::endl;
  }
  else msg_Tracking()<<METHOD<<" satisfied four-momentum conservation."<<std::endl;
#endif
  m_transitions += 1;
  return success;  
}


double Soft_Cluster_Handler::TransformWeight(Cluster * cluster,ATOOLS::Flavour & hadron,
					     const bool & lighter,const bool & enforce)
{
  Flavour_Pair fpair;
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;

  double MC(cluster->Mass());
  msg_Tracking()<<"@@ "<<METHOD<<" for cluster "<<cluster->Number()
		<<" ("<<fpair.first<<", "<<fpair.second<<", "<<"mass = "<<MC<<"): "<<endl
		<<"@@ Heavy = "<<p_singletransitions->GetHeaviestTransition(fpair)
		<<" ("<<p_singletransitions->GetHeaviestTransition(fpair).HadMass()<<"), "
		<<"light = "<<p_singletransitions->GetLightestTransition(fpair)
		<<" ("<<p_singletransitions->GetLightestTransition(fpair).HadMass()<<"), "
		<<"enforce = "<<enforce<<", lighter = "<<lighter<<std::endl;
  if (!enforce && p_singletransitions->GetHeaviestMass(fpair)<MC+m_transitionoffset) {
    msg_Tracking()<<"@@      --> too heavy, no transformation."<<std::endl;
    hadron = Flavour(kf_none);
    return 0.;
  }
  Single_Transition_Miter stiter = p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg_Tracking()<<"."<<std::endl;
    msg_Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<"\n"
	       <<"   will continue and hope for the best."<<endl;
    hadron = Flavour(kf_none);
    return 0.;
  }
  Single_Transition_List * stl(stiter->second);
  Single_Transition_Siter  start(stl->begin()),siter;
  if (lighter && hadron!=Flavour(kf_none)) {
    do {
      if (start->first==hadron) {
	siter = start;
	if ((++siter)!=stl->end()) start++;
	else return 0.;
	break;
      }
      else start++;
    } while (start!=stl->end());
  }

  double wt(0.),totweight(0.);
  
  for (siter=start;siter!=stl->end();siter++) {
    if (siter->first!=hadron) {
      wt  = TransformKin(MC,siter->first,enforce);
      wt *= siter->second;
    }
    totweight += wt;
  }

  double disc(totweight * 0.9999999999*ran.Get());
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

  msg_Tracking()<<"@@     --> "<<hadron<<" ("<<hadron.HadMass()<<", wt = "<<disc/totweight<<")."<<std::endl;

  return wt/(16.*M_PI*MC);
}

double Soft_Cluster_Handler::TransformKin(const double MC,const ATOOLS::Flavour & flav,
					  const bool & enforce) {
  double weight(1.);
  double mass(flav.HadMass()),mass2(mass*mass);
  double width(Max(flav.Width(),1.e-20)), width2(width*width);
  switch (m_tweightmode) {
  case TransWeight::waveonly:
    break;
  case TransWeight::nonrelativistic:
    if (!enforce && dabs(MC-mass)>20.*width) return 0.;
    weight *= pow(mass2/(sqr(MC-mass) + width2/4.),m_kappa);
    weight *= pow(width2/mass2,m_lambda);
    break;
  case TransWeight::relativistic:
  default:
    if (!enforce && dabs(MC-mass)>20.*width) return 0.;
    weight *= pow(sqr(mass2)/(sqr(MC*MC-mass2) + mass2*width2),m_kappa);
    weight *= pow(width2/mass2,m_lambda);
    break;
  }
  return weight;
}

double Soft_Cluster_Handler::DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2,
					 const bool enforce)
{
  if (m_dweightmode==DecayWeight::off && !enforce) return 0.;
  Flavour_Pair flpair;
  flpair.first  = cluster->GetTrip()->m_flav;
  flpair.second = cluster->GetAnti()->m_flav;

  Double_Transition_Miter dtliter = p_doubletransitions->GetTransitions()->find(flpair);
  double MC(cluster->Mass());
  msg_Tracking()<<"@@ "<<METHOD<<" for cluster "<<cluster->Number()
		<<" ("<<flpair.first<<", "<<flpair.second<<", "<<"mass = "<<MC<<"): "<<endl
		<<"@@ Heavy = "<<p_doubletransitions->GetHeaviestTransition(flpair).first
		<<" + "<<p_doubletransitions->GetHeaviestTransition(flpair).second
		<<" ("<<(p_doubletransitions->GetHeaviestTransition(flpair).first.HadMass()+
			 p_doubletransitions->GetHeaviestTransition(flpair).second.HadMass())<<"), "
		<<"light = "<<p_doubletransitions->GetLightestTransition(flpair).first
		<<" + "<<p_doubletransitions->GetLightestTransition(flpair).second
		<<" ("<<(p_doubletransitions->GetLightestTransition(flpair).first.HadMass()+
			 p_doubletransitions->GetLightestTransition(flpair).second.HadMass())<<"), "
		<<"enforce = "<<enforce<<"."<<std::endl;
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential Error in "<<METHOD<<" : "<<endl
	       <<"   No viable transition found for "<<flpair.first<<"/"<<flpair.second
	       <<" (mass = "<<MC<<") in "<<p_doubletransitions->GetTransitions()->size()<<"."<<endl
	       <<"   Return 'false' and hope for the best."<<std::endl;
    return 0.;
  }
  if (p_doubletransitions->GetLightestMass(flpair)>MC) {
    if (enforce) {
      msg_Tracking()<<"Warning in "<<METHOD<<" : "<<endl
		    <<"   No viable transition found for "<<flpair.first<<"/"<<flpair.second
		    <<" (mass = "<<MC<<")."<<endl
		    <<"   Return 'false' and hope for the best."<<std::endl;
    }
    return 0.;
  }
  if (!enforce && p_doubletransitions->GetHeaviestMass(flpair)<MC+m_decayoffset) {
    msg_Tracking()<<"@@      --> too heavy, no decay."<<std::endl;
    had1 = had2 = Flavour(kf_none);
    return 0.;
  }

  double totweight(0.),MC2(MC*MC),m1,m2,wt(1.);
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.HadMass();
    m2  = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = 1.;
      if (m_dweightmode==DecayWeight::phasespace ||
	  m_dweightmode==DecayWeight::phasespace_waves ||
	  m_dweightmode==DecayWeight::phasespace_masses_waves)
	wt *= sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
      if (m_dweightmode==DecayWeight::phasespace_masses_waves)
	wt *= pow((4.*m1*m2)/MC2,m_chi);
      if (m_dweightmode==DecayWeight::phasespace ||
	  m_dweightmode==DecayWeight::phasespace_waves ||
	  m_dweightmode==DecayWeight::waves)
	wt *= decit->second;
      totweight += wt;
    }
  }
  if (totweight<=0.) {
    msg_Tracking()<<"Error in "<<METHOD<<" :"<<std::endl
		  <<"   Cluster of mass "<<MC<<" from {"<<flpair.first<<", "
		  <<flpair.second<<"} passed mass conditions,"<<std::endl
		  <<"   but no viable transition found."<<std::endl
		  <<"   Return 0 and hope for the best."<<std::endl;
    return 0.;
  }

  had1 = had2 = Flavour(kf_none); 
  double disc(totweight * 0.9999999999*ran.Get());
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.HadMass();
    m2    = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = 1.;
      if (m_dweightmode==DecayWeight::phasespace ||
	  m_dweightmode==DecayWeight::phasespace_waves ||
	  m_dweightmode==DecayWeight::phasespace_masses_waves)
	wt *= sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
      if (m_dweightmode==DecayWeight::phasespace_masses_waves)
	wt *= pow((4.*m1*m2)/MC2,m_chi);
      if (m_dweightmode==DecayWeight::phasespace ||
	  m_dweightmode==DecayWeight::phasespace_waves ||
	  m_dweightmode==DecayWeight::waves)
	wt *= decit->second;
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
				      const Flavour had1,const Flavour had2)
{
  msg_Tracking()<<"## "<<METHOD<<" for "<<cluster->Number()<<" --> "<<had1<<" + "<<had2<<"."<<std::endl;
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.HadMass()), m22 = sqr(had2.HadMass());

  cluster->BoostInCMSAndRotateOnZ();
  int sign(cluster->GetTrip()->m_mom[3]<0?-1:1);

  double E1((M2+m12-m22)/(2.*M)), p2max(sqr(E1)-m12), pt2min(0.), pt(0.);
  if (cluster->GetTrip()->m_info=='L' || cluster->GetAnti()->m_info=='L') {
    pt2min       = Min(p2max*m_pt2maxfac/4.,m_pt2min/m12*m_pt2min/m22*m_pt2min);
    pt           = sqrt(p_as->SelectPT(p2max*m_pt2maxfac,pt2min));
  }
  else {
    // isotropic decay
    pt           = sqrt(p2max*m_pt2maxfac*(1.-sqr(1.-2.*ran.Get())));
  }
  //std::cout<<"Check this : pt = "<<pt<<" from p2_max =  "<<p2max<<" "
  //	   <<"and pt2_min = "<<pt2min<<"."<<std::endl<<(*cluster)<<std::endl;

  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,sign*pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

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

bool Soft_Cluster_Handler::ShiftMomenta(Cluster_List * clin)
{
  if (!TryLocalCompensation(clin)) return ForceMomenta(clin);
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
	  if ((cluster->GetNBTrip()->size()!=1 && cluster->GetNBAnti()->size()!=1) ||
	      (cluster->GetNBTrip()->size()==1 && cluster->GetNBAnti()->size()==1)) {
	    if (1.>0.5) {
	      partner = cluster->GetNBAnti();
	      direx   = -1;
	    }
	    else {
	      partner = cluster->GetNBTrip();
	      direx   = +1;
	    }
	  }
	  else if (cluster->GetNBTrip()->size()==1 && cluster->GetNBAnti()->size()!=1) {
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
	//msg_Out()<<"@@@ "<<METHOD<<": no partner for cluster "<<cluster->Number()<<"."<<std::endl;
	return false;
      }
      mass1 = (*cluster)[0].HadMass();
      mass2 = partner->size()==1?(*partner)[0].HadMass():partner->Mass();
      if (sqr(mass1+mass2)>(cluster->Momentum()+partner->Momentum()).Abs2()) {
	if (direx==0) {
	  //msg_Out()<<"@@@ "<<METHOD<<": Tried "<<partner->Number()
	  //	   <<" for cluster "<<cluster->Number()<<", not enough energy."<<std::endl;
	  return false;
	}
	if (direx==-1) partner = cluster->GetNBTrip();
	if (direx== 1) partner = cluster->GetNBAnti();
	mass2 = partner->size()==1?(*partner)[0].HadMass():partner->Mass();
	if (sqr(mass1+mass2)>(cluster->Momentum()+partner->Momentum()).Abs2()) {
	  //msg_Out()<<"@@@ "<<METHOD<<": Tried "<<partner->Number()
	  //	   <<" for cluster "<<cluster->Number()<<", not enough energy."<<std::endl;
	  return false;
	}
      }
      masses.clear();
      momenta.clear();
      masses.push_back(mass1);
      masses.push_back(mass2);
      momenta.push_back(cluster->Momentum());
      momenta.push_back(partner->Momentum());
      if (!hadpars.AdjustMomenta(2,&momenta.front(),&masses.front())) {
	//msg_Out()<<"@@@ "<<METHOD<<": Tried "<<partner->Number()
	//	 <<" for cluster "<<cluster->Number()<<", adjusting didn't work."<<std::endl;
	return false;
      }
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

  if (!hadpars.AdjustMomenta(size,&momenta.front(),&masses.front())) {
    if (size>1) {
      msg_Error()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : "<<std::endl
		 <<"   Could not adjust momenta for : "<<std::endl;
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Error()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Error()<<" ("<<((**cit)[0])<<" )";
	msg_Error()<<" for "<<std::endl<<(**cit)<<std::endl;
	i++;
      }
      msg_Error()<<"   Will possibly lead to retrying the event."<<std::endl;
    }
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
    msg_Error()<<METHOD<<" yields a momentum violation for  "<<size<<" : "<<std::endl
	       <<"   "<<checkbef<<" - "<<checkaft<<" --> "
	       <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl
	       <<(*clin)<<std::endl;
  }
  else msg_Tracking()<<METHOD<<" conserves momentum : "
		     <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl;
#endif
  return true;
}

Vec4D Soft_Cluster_Handler::SumMomentum(Cluster_List * clin) {
  Cluster_Iterator cit;
  Vec4D listmom(0.,0.,0.,0.);
  for (cit=clin->begin();cit!=clin->end();cit++) listmom += (*cit)->Momentum(); 
  return listmom;
}

