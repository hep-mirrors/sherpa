#include "Soft_Cluster_Handler.H"
#include "Random.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::Soft_Cluster_Handler() :
  m_stmode(stm::masswidthXwaves), 
  m_dtmode(dtm::waves_PS), 
  p_singletransitions(hadpars.GetSingleTransitions()),
  p_doubletransitions(hadpars.GetDoubleTransitions()),
  p_popper(hadpars.GetPopper()),
  m_offset1(hadpars.Get(string("Offset_C->H"))),
  m_offset2(hadpars.Get(string("Offset_C->HH"))),
  m_photonenergy(hadpars.Get(string("Photon_Energy")))
{ }

Soft_Cluster_Handler::~Soft_Cluster_Handler() {}

bool Soft_Cluster_Handler::TreatClusterList(Cluster_List * clin,Blob * blob)
{
  m_ctrans.clear();
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) 
    CheckCluster((*cit),clin->size()==1);

  m_ctit=m_ctrans.begin();
  if (clin->size()==1) {
    switch (m_ctit->second.size()) {
    case 1:
      msg.Error()<<"Potential problem in "<<METHOD<<endl;
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

  ShiftMomenta(clin->size());

  Flavour had1,had2;
  m_ctit=m_ctrans.begin();
  while (m_ctit!=m_ctrans.end()) {
    if (m_ctit->second.size()==1) {
      blob->AddToOutParticles(m_ctit->first->GetSelf());
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	if ((*cit)==m_ctit->first) { 
	  delete (*cit);
	  clin->erase(cit);
	  break;
	}
      }
    }
    else if (m_ctit->second.size()==2) {
      FixHHDecay(m_ctit->first,m_ctit->second[0],m_ctit->second[1]);
    }
    m_ctrans.erase(m_ctit++);
  }
  return true;
}

void Soft_Cluster_Handler::CheckCluster(Cluster * cluster,bool mustdecay)
{
  //   cout<<"----------------------------------------------------------------------"<<endl
  //       <<"Check this out : "<<cluster->Mass()<<"("<<mustdecay<<") / "
  //       <<cluster->GetFlav(1)<<" + "<<cluster->GetFlav(2)<<endl;
  vector<Flavour> flist;
  m_ctrans[cluster] = flist;
  Flavour had1,had2,hadron;
  double decayweight     = DecayWeight(cluster,had1,had2); 
  double transformweight = TransformWeight(cluster,hadron); 
  if (decayweight<=0. && transformweight<=0.) return;

  // single cluster in singlet list
  if (mustdecay==true) {
    if (decayweight/(decayweight+transformweight)>ran.Get()) {
      // decay is favoured -> no problem
      m_ctrans[cluster].push_back(had1);
      m_ctrans[cluster].push_back(had2);
    }
    else {
      if (hadron.Mass()+m_offset2<cluster->Mass()) {
	// selected hadron for transformation is lighter, compensate with photon
	m_ctrans[cluster].push_back(hadron);
	m_ctrans[cluster].push_back(Flavour(kf::photon));
      }
      else {
	// selected hadron for transformation is heavier, select a new one
	transformweight = TransformWeight(cluster,hadron,true); 
	if (decayweight/(decayweight+transformweight)>ran.Get()) {
	  // but now decay is favoured -> no problem
	  m_ctrans[cluster].push_back(had1);
	  m_ctrans[cluster].push_back(had2);
	}
	else {
	  // force decay into lighter hadron plus photon, irrespecitve of how soft.
	  m_ctrans[cluster].push_back(hadron);
	  m_ctrans[cluster].push_back(Flavour(kf::photon));
	}
      }
    }
    //     cout<<METHOD<<" ("<<mustdecay<<", "<<m_ctrans[cluster].size()<<") : "
    //      	<<decayweight<<" ("<<had1<<","<<had2<<") ; "
    //      	<<transformweight<<" ("<<hadron<<") for "
    //      	<<cluster->Mass()<<" / "<<cluster->GetFlav(1)<<" + "<<cluster->GetFlav(2)<<endl
    //      	<<"----------------------------------------------------------------------"<<endl;
    return;
  }
  // regular case.
  if (decayweight/(decayweight+transformweight)>ran.Get()) {
    m_ctrans[cluster].push_back(had1);
    m_ctrans[cluster].push_back(had2);
  }
  else m_ctrans[cluster].push_back(hadron);

  //   cout<<METHOD<<" ("<<mustdecay<<") : "<<decayweight<<" ("<<had1<<","<<had2<<") ; "
  //       <<transformweight<<" ("<<hadron<<") for "
  //       <<cluster->Mass()<<" / "<<cluster->GetFlav(1)<<" + "<<cluster->GetFlav(2)<<endl
  //       <<"----------------------------------------------------------------------"<<endl;
  return;
}
 

double Soft_Cluster_Handler::DecayWeight(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  had1 = had2 = Flavour(kf::none);
  Flavour_Pair flpair;
  flpair.first  = cluster->GetFlav(1);
  flpair.second = cluster->GetFlav(2);

  Double_Transition_Miter dtliter = p_doubletransitions->GetTransitions()->find(flpair);
  if (dtliter==p_doubletransitions->GetTransitions()->end()) {
    msg.Error()<<"ERROR in "<<METHOD<<" : "<<endl
	       <<"   No transition table found for "<<flpair.first<<"/"<<flpair.second<<endl
	       <<"   Return 'false' and hope for the best."<<std::endl;
    return -1.;
  }

  double totweight(0.), MC(cluster->Mass()), MC2(MC*MC);
  double wt,m1,m2;
  //   cout<<METHOD<<" : "<<flpair.first<<"/"<<flpair.second<<" -> "
  //       <<p_doubletransitions->GetHeaviestMass(flpair)<<","
  //       <<p_doubletransitions->GetLightestMass(flpair)<<" vs. "<<MC<<endl;

  if (p_doubletransitions->GetHeaviestMass(flpair)+m_offset2<MC) return 0.;
  if (p_doubletransitions->GetLightestMass(flpair)>MC)           return 0.;
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.Mass();
    m2  = decit->first.second.Mass();
    if (m1+m2<MC) {
      wt = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
      if (m_dtmode==dtm::waves_PS) wt *= decit->second;
      totweight += wt;
    }
  }
  if (totweight==0.) return 0.;

  had1 = had2 = Flavour(kf::none); 
  double disc = totweight*ran.Get();
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.Mass();
    m2    = decit->first.second.Mass();
    if (m1+m2>MC) continue;
    wt = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2)));
    if (m_dtmode==dtm::waves_PS) wt *= decit->second;
    disc -= wt;
    if (disc<0.) {
      had1 = decit->first.first;
      had2 = decit->first.second;
      break;
    }
  }
  //  cout<<METHOD<<":"<<totweight/(16.*M_PI*MC*MC*MC)<<" --> "<<had1<<" & "<<had2<<endl;
  return totweight * 1./(16.*M_PI*MC*MC*MC);
}

void Soft_Cluster_Handler::FixHHDecay(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.PSMass()), m22 = sqr(had2.PSMass());
  double ptmax   = sqrt(sqr(M2-m12-m22)-4.*m12*m22)/(2.*M); 

  double pt      = p_popper->SelectPT(ptmax);

  cluster->BoostInCMSAndRotateOnZ();
  double E1      = (M2+m12-m22)/(2.*M);
  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  Particle * part;
  Cluster * clus;
  clus = new Cluster();
  clus->SetMomentum(0,p1);
  clus->SetPrev(cluster);
  cluster->SetLeft(clus);

  clus = new Cluster();
  clus->SetMomentum(0,p2);
  clus->SetPrev(cluster);
  cluster->SetRight(clus);

  cluster->RotateAndBoostBack();

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
  fpair.first  = cluster->GetFlav(1);
  fpair.second = cluster->GetFlav(2);
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return -1.;
  Single_Transition_Miter stiter = p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg.Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<";"<<endl<<(*cluster->GetPrev())
	       <<"   will continue and hope for the best."<<endl;
    return -1.;
  }

  double totweight(0.), MC(cluster->Mass()), MC2(MC*MC);
  double wt,m1,m2;
  //cout<<METHOD<<" : "<<fpair.first<<"/"<<fpair.second<<" -> "
  //   <<p_singletransitions->GetHeaviestMass(fpair)<<","
  //   <<p_singletransitions->GetLightestMass(fpair)<<" vs. "<<MC<<endl;
  if (p_singletransitions->GetHeaviestMass(fpair)+m_offset1<MC) return 0.;
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
  Hadron_Wave_Function * waves;
  Single_Transition_Siter start=stl->begin();
  if (lighter && hadron!=Flavour(kf::none)) {
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
  if (lighter && hadron!=Flavour(kf::none)) {
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
  // cout<<METHOD<<" --> "<<totweight/(2.*MC)<<" for "<<hadron<<endl;
  return totweight/(2.*MC);;
}

void Soft_Cluster_Handler::ShiftMomenta(const int size)
{
  double masses[size];
  Vec4D  momenta[size];
  int    pos(0);
  for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
    if (m_ctit->second.size()==1) {
      masses[pos]  = m_ctit->second.begin()->Mass();
      momenta[pos] = m_ctit->first->Momentum();
      pos++;
    }
  }
  if (pos<2) {
    for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
      if (m_ctit->second.size()!=1) {
	masses[pos]  = m_ctit->first->Mass();
	momenta[pos] = m_ctit->first->Momentum();
	pos++;
      }
    }
  }
  hadpars.AdjustMomenta(pos,momenta,masses);
  pos = 0;
  for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
    if (m_ctit->second.size()==1) {
      m_ctit->first->SetMomentum(0,momenta[pos]);
      m_ctit->first->GetSelf()->SetMomentum(momenta[pos]);
      m_ctit->first->GetSelf()->SetFinalMass(masses[pos]);
      m_ctit->first->GetSelf()->SetFlav((*m_ctit->second.begin()));
      pos++;
    }
  }
  if (pos<2) {
    for (m_ctit=m_ctrans.begin();m_ctit!=m_ctrans.end();m_ctit++) {
      if (m_ctit->second.size()!=1) {
	m_ctit->first->RescaleMomentum(momenta[pos]);
	pos++;
      }
    }
  }
}
