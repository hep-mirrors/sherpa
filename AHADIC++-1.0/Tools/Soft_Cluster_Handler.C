#include "Soft_Cluster_Handler.H"
#include "Random.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler(Single_Transitions * singletransitions,
					   Double_Transitions * doubletransitions,
					   const double offset1,const double offset2,
					   const double alpha,const double eta,
					   const double photonenergy) :
  m_stmode(stm::masswidthXwaves), m_dtmode(dtm::waves_PS), 
  p_singletransitions(singletransitions), p_doubletransitions(doubletransitions),
  m_offset1(offset1), m_offset2(offset2),
  m_alpha(alpha), m_eta(eta), m_photonenergy(photonenergy)
{ }

Soft_Cluster_Handler::~Soft_Cluster_Handler() {}

bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin, Blob * blob)
{
  Cluster_Iterator cit;
  Cluster * cluster;

  int size(0);
  for (cit=clin->begin();cit!=clin->end();cit++) {
    if ((*cit)==NULL || !(*cit)->Active()) continue;
    size += CheckCluster((*cit),clin->size()==1);
  }
  if (size==0) return true;

  if (clin->size()==1) {
    cluster = (*clin->begin());
    switch (cluster->size()) {
    case 1:
      msg_Error()<<"Potential problem in "<<METHOD<<" : "<<std::endl
		 <<"   Clusterlist with one element that needs to transform to a hadron."<<std::endl
		 <<"   Will possibly lead to retrying the event."<<std::endl;
      return false;
    case 2:
      FixHHDecay(cluster,blob);
      cluster->SetActive(false);
      clin->clear();
    case 0:
    default:
      return true;
    }
  }

#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) checkbef += (*cit)->Momentum(); 
#endif

  double E(-1.);
  bool breakit(true);
  while (!CheckIfAllowed(clin,clin->size(),E)) {
    breakit = UpdateTransitions(clin,clin->size());
    if (!breakit) break;
  }

  if (!ShiftMomenta(clin,clin->size(),!breakit)) {
      msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
		 <<"   Could not shift momenta."<<std::endl
		 <<"   Will possibly lead to retrying the event."<<std::endl;
    return false;
  }

  Flavour had1,had2;
  for (cit=clin->begin();cit!=clin->end();) {
    cluster = (*cit);
    switch (cluster->size()) {
    case 2:
      FixHHDecay(cluster,blob);
#ifdef AHAmomcheck
      checkaft += cluster->GetLeft()->Momentum();
      checkaft += cluster->GetRight()->Momentum();
#endif
      cluster->SetActive(false);
      cit = clin->erase(cit);
      break;
    case 1:
      blob->AddToOutParticles(cluster->GetSelf());
#ifdef AHAmomcheck
      checkaft += cluster->Momentum();
      msg_Debugging()<<METHOD<<" involving a C->H transition."<<std::endl;
#endif
      cluster->SetActive(false);
      cit = clin->erase(cit);
      break;
    case 0:
#ifdef AHAmomcheck
      checkaft += cluster->Momentum();
#endif
      cit++;
      break;
    }
  }

#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif
  return true;
}

int Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool mustdecay,const double Mmax)
{
  cluster->clear();
  Flavour had1,had2,hadron;

  double decayweight(0.), transformweight(0.);
  if (p_singletransitions->MustDesintegrate(cluster,had1,had2)) {
    decayweight=1.;
  }
  else {
    decayweight     = DecayWeight(cluster,had1,had2);
    transformweight = TransformWeight(cluster,hadron); 
  }
  if (decayweight<=0. && transformweight<=0.) {
    return 0;
  }

  // single cluster in singlet list
  if (mustdecay==true) {
    if (decayweight/(decayweight+transformweight)>ran.Get()) {
      cluster->push_back(had1);
      cluster->push_back(had2);
    }
    else {
      if (hadron.Mass()+m_photonenergy<cluster->Mass()) {
	// selected hadron for transformation is lighter, compensate with photon
	cluster->push_back(hadron);
	cluster->push_back(Flavour(kf_photon));
      }
      else {
	// selected hadron for transformation is heavier, select a new one
	transformweight = TransformWeight(cluster,hadron,true); 
	if (decayweight/(decayweight+transformweight)>ran.Get()) {
	  // but now decay is favoured -> no problem
	  cluster->push_back(had1);
	  cluster->push_back(had2);
	}
	else {
	  // force decay into lighter hadron plus photon, irrespecitve of how soft.
	  cluster->push_back(hadron);
	  cluster->push_back(Flavour(kf_photon));
	}
      }
    }
    return 2;
  }
  // regular case.
  if (decayweight/(decayweight+transformweight)>ran.Get()) {
    cluster->push_back(had1);
    cluster->push_back(had2);
    return 2;
  }
  else cluster->push_back(hadron);
  return 1;
}
 
double Soft_Cluster_Handler::DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2,
					 const double Mmax)
{
  had1 = had2 = Flavour(kf_none);
  Flavour_Pair flpair;
  flpair.first  = cluster->GetTrip()->m_flav;
  flpair.second = cluster->GetAnti()->m_flav;

  Double_Transition_Miter dtliter = p_doubletransitions->GetTransitions()->find(flpair);
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg_Error()<<"ERROR in "<<METHOD<<" : "<<endl
	       <<"   No transition table found for "<<flpair.first<<"/"<<flpair.second<<endl
	       <<"   Return 'false' and hope for the best."<<std::endl;
    return 0.;
  }

  double totweight(0.), MC(cluster->Mass()), MC2(MC*MC);
  double wt,m1,m2;
  bool   forceit(false);
  double altweight(0.);
  if (flpair.first.Mass()+flpair.second.Mass()+
      2.*p_singletransitions->GetLightestConstituent().Mass()>MC) {
    forceit=true;
  }
  else {
    if (p_doubletransitions->GetHeaviestMass(flpair)<MC+m_offset2) return 0.;
    if (p_doubletransitions->GetLightestMass(flpair)>MC)           return 0.;
  }
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.Mass();
    m2  = decit->first.second.Mass();
    if (m1+m2<MC) {
      wt = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * pow(4.*m1*m2/MC2,m_alpha);
      if (m_dtmode==dtm::waves_PS) wt *= decit->second;
      altweight += wt;
      if  (m1+m2<MC+int(forceit)*m_offset2) totweight += wt;
    }
  }
  if (totweight==0. && altweight==0.) return 0.;

  had1 = had2 = Flavour(kf_none); 
  double disc = totweight*ran.Get();
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.Mass();
    m2    = decit->first.second.Mass();
    if (m1+m2<MC && m1+m2<MC+int(forceit)*m_offset2) {
      wt = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * pow(4.*m1*m2/MC2,m_alpha);
      if (m_dtmode==dtm::waves_PS) wt *= decit->second;
      disc -= wt;
      if (disc<0.) {
	had1 = decit->first.first;
	had2 = decit->first.second;
	break;
      }
    }
  }

  Flavour alt1 = Flavour(kf_none), alt2 = Flavour(kf_none);
  if (altweight!=0) {
    disc = altweight*ran.Get();
    for (Double_Transition_Siter decit=dtliter->second->begin();
	 decit!=dtliter->second->end();decit++) {
      m1    = decit->first.first.Mass();
      m2    = decit->first.second.Mass();
      if (m1+m2<MC) {
	wt = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * pow(4.*m1*m2/MC2,m_alpha);
	if (m_dtmode==dtm::waves_PS) wt *= decit->second;
	disc -= wt;
	if (disc<0.) {
	  alt1 = decit->first.first;
	  alt2 = decit->first.second;
	  break;
	}
      }
    }
  }
  return totweight * 1./(16.*M_PI*MC*MC*MC);
}

void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Blob * blob)
{
#ifdef AHAmomcheck
  Vec4D  checkbef = cluster->Momentum();
#endif
  Flavour had1((*cluster)[0]), had2((*cluster)[1]);

  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.PSMass()), m22 = sqr(had2.PSMass());

  cluster->BoostInCMSAndRotateOnZ();
  double E1      = (M2+m12-m22)/(2.*M);
  double sinthet = pow(ran.Get(),m_eta);
  double pt      = sqrt(sqr(E1)-m12)*sinthet;
  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  Particle * part;
  Cluster * clus;
  clus = new Cluster(p1,had1,false);
  clus->SetPrev(cluster);
  cluster->SetLeft(clus);
#ifdef memchecker
  std::cout<<"@@@ New cluster "<<clus<<" from "<<METHOD<<"."<<std::endl;
#endif

  clus = new Cluster(p2,had2,false);
  clus->SetPrev(cluster);
  cluster->SetRight(clus);
#ifdef memchecker
  std::cout<<"@@@ New cluster "<<clus<<" from "<<METHOD<<"."<<std::endl;
#endif

  cluster->RotateAndBoostBack();

  blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
  blob->AddToOutParticles(cluster->GetRight()->GetSelf());


#ifdef AHAmomcheck
  if (dabs((checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2())>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"    Four-momentum not conserved: "
	       <<checkbef<<" vs. "<<(cluster->GetLeft()->Momentum()+cluster->GetRight()->Momentum())
	       <<" : "
	       <<(checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2()
	       <<"."<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum."<<std::endl;
#endif

}






double Soft_Cluster_Handler::TransformWeight(Cluster * cluster,ATOOLS::Flavour & hadron,
					     const double Mmax,bool lighter)
{
  Flavour_Pair fpair;
  double tfwt(0.);
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;
  Single_Transition_Miter stiter = p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<";"<<endl<<(*cluster->GetPrev())
	       <<"   will continue and hope for the best."<<endl;
    return tfwt;
  }

  double MC(Mmax<0?cluster->Mass():Mmax);
  if (p_doubletransitions->GetLightestMass(fpair)<MC) {
    if (p_singletransitions->GetHeaviestMass(fpair)+m_offset1<MC && !lighter) {
      return tfwt;
    }
  }
  
  switch (m_stmode) {
  case (stm::massXwaves) :
  case (stm::simplemass) :
    tfwt = SimpleMassCriterion(stiter->second,hadron,MC,lighter);
    break;
  case (stm::masswidthXwaves) :
  case (stm::masswidth) :
  default :
    tfwt = MWCriterion(stiter->second,hadron,MC,lighter);
  }
  return tfwt;
}

double Soft_Cluster_Handler::SimpleMassCriterion(Single_Transition_List * stl,
						 Flavour & hadron,
						 const double MC, bool lighter)
{
  double wt(0.),maxwt(0.);
  Single_Transition_Siter start=stl->begin();
  if (lighter && hadron!=Flavour(kf_none)) {
    for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
      if (siter->first==hadron) {
	start = siter; 
	if ((siter++)!=stl->end()) start++;
	if (dabs(start->first.Mass()-hadron.Mass())<1.e-2) {
	  if ((siter++)!=stl->end()) start++;
	}
	break;
      }
    }
  }
  

  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (lighter && siter->first.Mass()>MC) continue;
    wt = 1./pow(sqr(sqr(MC)-sqr(siter->first.Mass()))+1.e-8,0.25);
    if (m_stmode==stm::massXwaves) wt *= siter->second;
    if (wt>maxwt) {
      hadron = siter->first;
      maxwt  = wt;
    }
  }
  return maxwt;
}

double Soft_Cluster_Handler::MWCriterion(Single_Transition_List * stl,
					 Flavour & hadron,
					 const double MC, bool lighter)
{
  double wt(0.), totweight(0.), minwidth(1.e6);
  Single_Transition_Siter start=stl->begin();
  if (lighter && hadron!=Flavour(kf_none)) {
    for (Single_Transition_Siter siter=stl->begin();siter!=stl->end();siter++) {
      if (siter->first==hadron) {
	start = siter;
	if ((siter++)!=stl->end()) start++;
	if (start->first.Mass()-hadron.Mass()<1.e-2) {
	  if ((siter++)!=stl->end()) start++;
	}
	break;
      }
    }
  }
  
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (!siter->first.IsStable() && (siter->first.Width()<minwidth)) 
      minwidth = siter->first.Width();
  }
  if (minwidth==1.e6) return SimpleMassCriterion(stl,hadron,MC,lighter);
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (lighter && siter->first.Mass()>MC) continue;
    if (siter->first.IsStable()) 
      wt = sqr(siter->first.Mass()*minwidth)/
	(sqr(sqr(MC)-sqr(siter->first.Mass())) + sqr(siter->first.Mass()*minwidth));
    else 
      wt = sqr(siter->first.Mass()*siter->first.Width())/
	(sqr(sqr(MC)-sqr(siter->first.Mass())) + sqr(siter->first.Mass()*siter->first.Width()));
    if (m_stmode==stm::masswidthXwaves) wt *= siter->second;
    if (wt>totweight) {
      hadron    = siter->first;
      totweight = wt;
    }
  }
  return totweight/(2.*MC);;
}

bool Soft_Cluster_Handler::CheckIfAllowed(Cluster_List * clin,int size,double & E) {
  double totmass(0.);
  Vec4D  totmom(0.,0.,0.,0.);
  Cluster * cluster;
  int count(0), dec(0);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      totmass += (*cluster)[0].Mass();
      if (E<0) totmom += cluster->Momentum();
      count++;
    }
    else if (cluster->size()==2) dec++;
  }
  if (count<1) return true;
  if (E<0) E = sqrt(totmom.Abs2());
  return (totmass<E);
}

bool Soft_Cluster_Handler::UpdateTransitions(Cluster_List * clin,int size) {
  Cluster * cluster;
  Flavour hadron;
  bool    transform(false);
  double  tfwt;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      hadron = (*cluster)[0]; 
      tfwt   = TransformWeight(cluster,hadron,hadron.Mass(),true);
      if (tfwt>0.) {
	cluster->clear();
	cluster->push_back(hadron);
	transform = true;
      }
    }
  }
  return transform;
}

bool Soft_Cluster_Handler::ShiftMomenta(Cluster_List * clin,int size,bool takeall)
{
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.),checkaft(0.,0.,0.,0.);
#endif
  std::vector<double> masses(size);
  std::vector<Vec4D>  momenta(size);
  int    pos(0),pos1;
  Cluster * cluster;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      masses[pos]  = (*cluster)[0].Mass();
      momenta[pos] = cluster->Momentum();
#ifdef AHAmomcheck
      checkbef    += momenta[pos];
#endif
      pos++;
    }
  }
  if (pos<2 || takeall) {
    for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
      cluster = (*cit);
      if (cluster->size()!=1) {
	masses[pos]  = cluster->Mass();
	momenta[pos] = cluster->Momentum();
#ifdef AHAmomcheck
	checkbef    += momenta[pos];
#endif 
	pos++;
      }
    }
  }
  
  if (!hadpars.AdjustMomenta(pos,&momenta.front(),&masses.front())) {
    if (pos>1) {
      msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
		 <<"   Could not adjust momenta for : "<<std::endl;
      for (int i=0;i<pos;i++) {
	msg_Error()<<"   "<<i<<" th particle: "<<masses[i]<<" "<<momenta[i]<<std::endl;
      }
      msg_Error()<<"   Will possibly lead to retrying the event."<<std::endl;
    }
    return false;
  }
  pos = 0;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
#ifdef AHAmomcheck
      checkaft += momenta[pos];
#endif
      cluster->SetMomentum(momenta[pos]);
      cluster->SetFlav((*cluster)[0]);
      pos++;
    }
  }
  pos1 = pos;
  if (pos<2 || takeall) {
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
      if (cluster->size()!=1) {
	cluster->RescaleMomentum(momenta[pos]);
#ifdef AHAmomcheck
	checkaft += cluster->Momentum();
#endif
	pos++;
      }
    }
  }
#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields a momentum violation for  "<<pos1<<" : "
  	     <<checkbef<<" - "<<checkaft<<" --> "
  	     <<(checkbef-checkaft).Abs2()<<"("<<pos1<<", "<<size<<")."<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum : "
		      <<(checkbef-checkaft).Abs2()<<"("<<pos1<<", "<<size<<")."<<std::endl;
#endif
  return true;
}
