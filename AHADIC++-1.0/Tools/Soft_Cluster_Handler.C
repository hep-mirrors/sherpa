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

bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin,Blob * blob)
{
  m_ctrans.clear();
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    if ((*cit)==NULL) continue;
    CheckCluster((*cit),clin->size()==1);
  }


  m_ctit=m_ctrans.begin();
  PRINT_VAR(clin->size());
  PRINT_VAR(*(m_ctit->first));
  if (clin->size()==1) {
    PRINT_VAR(m_ctit->second.size());
    PRINT_VAR(*(m_ctit->first));
    switch (m_ctit->second.size()) {
    case 1:
      msg_Error()<<"Potential problem in "<<METHOD<<endl;
      m_ctrans.clear();
      return false;
    case 2:
      FixHHDecay(m_ctit->first,m_ctit->second[0],m_ctit->second[1]);
    case 0:
    default:
      m_ctrans.clear();
      return true;
    }
  }
  PRINT_VAR(*(m_ctit->first));

  ShiftMomenta(clin->size());

#ifdef AHAmomcheck
    Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
    for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) checkbef += (*cit)->Momentum(); 
#endif
  Flavour had1,had2;
  m_ctit=m_ctrans.begin();
  while (m_ctrans.size()>0) {
    if (m_ctit->second.size()==1) {
      blob->AddToOutParticles(m_ctit->first->GetSelf());
#ifdef AHAmomcheck
      checkaft += m_ctit->first->Momentum();
      if (dabs((m_ctit->first->GetSelf()->Momentum()-m_ctit->first->Momentum()).Abs2())>1.e-12) {
	msg_Out()<<METHOD<<" cluster and particle momentum do not agree: "
		 <<(m_ctit->first->GetSelf()->Momentum()-m_ctit->first->Momentum())
		 <<" ("<<(m_ctit->first->GetSelf()->Momentum()-m_ctit->first->Momentum()).Abs2()<<")."
		 <<std::endl;
      }
#endif
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	if ((*cit)==m_ctit->first) { 
	  delete (*cit);
	  clin->erase(cit);
#ifdef AHAmomcheck
	  msg_Out()<<METHOD<<" involving a C->H transition."<<std::endl;
#endif
	  break;
	}
      }
    }
    else if (m_ctit->second.size()==2) {
      FixHHDecay(m_ctit->first,m_ctit->second[0],m_ctit->second[1]);
#ifdef AHAmomcheck
      checkaft += m_ctit->first->GetLeft()->Momentum();
      checkaft += m_ctit->first->GetRight()->Momentum();
#endif
    }
    if (m_ctit==m_ctrans.end()) break; 
    m_ctrans.erase(m_ctit++);
  }

#ifdef AHAmomcheck
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    if (!(*cit)->GetLeft()) checkaft += (*cit)->Momentum(); 
  }
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields momentum violation : "<<std::endl
	     <<checkbef<<" - "<<checkaft<<" --> "<<(checkbef-checkaft).Abs2()<<std::endl;
  }
  else msg_Out()<<METHOD<<" conserves momentum."<<std::endl;
#endif
  return true;
}

void Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool mustdecay)
{
  vector<Flavour> flist;
  if(m_ctrans.begin()!=m_ctrans.end()) PRINT_VAR(*m_ctrans.begin()->first);
  PRINT_VAR(*cluster);
  m_ctrans[cluster] = flist;
  Flavour had1,had2,hadron;

  double decayweight(0.), transformweight(0.);
  PRINT_VAR(*m_ctrans.begin()->first);
  PRINT_VAR((m_ctrans.begin()->first==cluster));
  if (p_singletransitions->MustDesintegrate(cluster,had1,had2)) {
    decayweight=1.;
    PRINT_VAR(*m_ctrans.begin()->first);
  }
  else {
    decayweight     = DecayWeight(cluster,had1,had2);
    PRINT_VAR(*m_ctrans.begin()->first);
    transformweight = TransformWeight(cluster,hadron); 
    PRINT_VAR(*m_ctrans.begin()->first);
  }
  if (decayweight<=0. && transformweight<=0.) {
    return;
  }

  // single cluster in singlet list
  if (mustdecay==true) {
    if (decayweight/(decayweight+transformweight)>ran.Get()) {
      // decay is favoured -> no problem
      m_ctrans[cluster].push_back(had1);
      m_ctrans[cluster].push_back(had2);
    }
    else {
      if (hadron.Mass()+m_photonenergy<cluster->Mass()) {
	// selected hadron for transformation is lighter, compensate with photon
	m_ctrans[cluster].push_back(hadron);
	m_ctrans[cluster].push_back(Flavour(kf_photon));
      }
      else {
	// selected hadron for transformation is heavier, select a new one
	transformweight = TransformWeight(cluster,hadron,true); 
        PRINT_VAR(*m_ctrans.begin()->first);
	if (decayweight/(decayweight+transformweight)>ran.Get()) {
	  // but now decay is favoured -> no problem
	  m_ctrans[cluster].push_back(had1);
	  m_ctrans[cluster].push_back(had2);
	}
	else {
	  // force decay into lighter hadron plus photon, irrespecitve of how soft.
	  m_ctrans[cluster].push_back(hadron);
	  m_ctrans[cluster].push_back(Flavour(kf_photon));
	}
      }
    }
    return;
  }
  // regular case.
  if (decayweight/(decayweight+transformweight)>ran.Get()) {
    m_ctrans[cluster].push_back(had1);
    m_ctrans[cluster].push_back(had2);
  }
  else m_ctrans[cluster].push_back(hadron);
  PRINT_VAR(*m_ctrans.begin()->first);

  return;
}
 
double Soft_Cluster_Handler::DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  PRINT_VAR(*m_ctrans.begin()->first);
  PRINT_VAR(*cluster);
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
  PRINT_VAR(*m_ctrans.begin()->first);

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
  PRINT_VAR(*m_ctrans.begin()->first);
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
  PRINT_VAR(*m_ctrans.begin()->first);

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
  PRINT_VAR(*m_ctrans.begin()->first);
  return totweight * 1./(16.*M_PI*MC*MC*MC);
}

void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Flavour & had1,Flavour & had2)
{
#ifdef AHAmomcheck
  Vec4D checkbef = cluster->Momentum();
#endif

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
  clus = new Cluster();
  clus->SetMomentum(p1);
  clus->SetPrev(cluster);
  cluster->SetLeft(clus);

  clus = new Cluster();
  clus->SetMomentum(p2);
  clus->SetPrev(cluster);
  cluster->SetRight(clus);

  cluster->RotateAndBoostBack();

#ifdef AHAmomcheck
  if (dabs((checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2())>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"    Four-momentum not conserved: "
	       <<checkbef<<" vs. "<<(cluster->GetLeft()->Momentum()+cluster->GetRight()->Momentum())
	       <<" : "<<(checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2()
	       <<"."<<std::endl;
  }
  else msg_Out()<<METHOD<<" conserves momentum."<<std::endl;
#endif

  part = new Particle(-1,had1,cluster->GetLeft()->Momentum());
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had1.PSMass());
  control::s_AHAparticles++;
  cluster->GetLeft()->SetSelf(part);
  
  part = new Particle(-1,had2,cluster->GetRight()->Momentum());
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had2.PSMass());
  control::s_AHAparticles++;
  cluster->GetRight()->SetSelf(part);
}

double Soft_Cluster_Handler::TransformWeight(Cluster * cluster,ATOOLS::Flavour & hadron,
					     bool lighter)
{
  Flavour_Pair fpair;
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;
  Single_Transition_Miter stiter = p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<";"<<endl<<(*cluster->GetPrev())
	       <<"   will continue and hope for the best."<<endl;
    return 0.;
  }

  double MC(cluster->Mass());
  if (p_doubletransitions->GetLightestMass(fpair)<MC) {
    if (p_singletransitions->GetHeaviestMass(fpair)+m_offset1<MC) return 0.;
  }
  switch (m_stmode) {
  case (stm::massXwaves) :
  case (stm::simplemass) :
    return SimpleMassCriterion(stiter->second,hadron,MC,lighter);
    break;
  case (stm::masswidthXwaves) :
  case (stm::masswidth) :
  default :
    return MWCriterion(stiter->second,hadron,MC,lighter);
  }
  return 0.;
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
	if (start->first.Mass()-hadron.Mass()<1.e-2) {
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
  
  //  cout<<METHOD<<" start with "<<start->first<<" / "<<hadron<<endl;
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
    //cout<<METHOD<<" "<<wt<<" "<<siter->second<<endl;
    if (wt>totweight) {
      hadron    = siter->first;
      totweight = wt;
    }
  }
  //cout<<METHOD<<" --> "<<totweight/(2.*MC)<<" for "<<hadron<<endl;
  return totweight/(2.*MC);;
}

void Soft_Cluster_Handler::ShiftMomenta(const int size)
{
  std::vector<double> masses(size);
  std::vector<Vec4D>  momenta(size);
  int    pos(0),pos1;
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.),checkaft(0.,0.,0.,0.);
#endif
  for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
    if (m_ctit->second.size()==1) {
      masses[pos]  = m_ctit->second.begin()->Mass();
      momenta[pos] = m_ctit->first->Momentum();
      PRINT_VAR(pos<<" "<<momenta[pos]);
#ifdef AHAmomcheck
      checkbef    += momenta[pos];
#endif
      pos++;
    }
  }
  if (pos<2) {
    for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
      if (m_ctit->second.size()!=1) {
	masses[pos]  = m_ctit->first->Mass();
	momenta[pos] = m_ctit->first->Momentum();
        PRINT_VAR(*(m_ctit->first));
        PRINT_VAR(pos<<" "<<momenta[pos]);
#ifdef AHAmomcheck
	checkbef    += momenta[pos];
#endif 
	pos++;
      }
    }
  }

  PRINT_VAR(pos<<" "<<size<<" "<<momenta.size());
  hadpars.AdjustMomenta(pos,&momenta.front(),&masses.front());

  pos = 0;
  for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
    if (m_ctit->second.size()==1) {
#ifdef AHAmomcheck
      checkaft += momenta[pos];
#endif
      m_ctit->first->SetMomentum(momenta[pos]);
      m_ctit->first->GetSelf()->SetMomentum(momenta[pos]);
      m_ctit->first->GetSelf()->SetFinalMass(masses[pos]);
      m_ctit->first->GetSelf()->SetFlav((*m_ctit->second.begin()));
      pos++;
    }
  }
  pos1 = pos;
  if (pos<2) {
    for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
      if (m_ctit->second.size()!=1) {
        PRINT_VAR(pos<<" "<<momenta.size()<<" "<<momenta[pos]);
	m_ctit->first->RescaleMomentum(momenta[pos]);
#ifdef AHAmomcheck
	checkaft += m_ctit->first->Momentum();
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
  else msg_Out()<<METHOD<<" conserves momentum : "
  		<<(checkbef-checkaft).Abs2()<<"("<<pos1<<", "<<size<<")."<<std::endl;
#endif
}
