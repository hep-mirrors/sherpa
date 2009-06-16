#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Soft_Cluster_Handler::
Soft_Cluster_Handler(Strong_Coupling * as,
		     Single_Transitions * singletransitions,
		     Double_Transitions * doubletransitions,
		     const double offset1,const double offset2,
		     const double kappa,
		     const double alpha,const double eta,
		     const double photonenergy,bool ana) :
  p_as(as),
  p_singletransitions(singletransitions), 
  p_doubletransitions(doubletransitions),
  m_ptmode(PTdist::alphadist),m_offset1(offset1), m_offset2(offset2),
  m_kappa(kappa), m_alpha(alpha), m_eta(eta), m_photonenergy(photonenergy),
  m_ana(ana)
{ 
  //std::cout<<METHOD<<":"<<m_kappa<<std::endl;exit(1);
  if (m_ana) {
    m_histograms[string("PT_HH")]  = new Histogram(0,0.,2.0,100);
    m_histograms[string("PT2_HH")] = new Histogram(0,0.,4.0,200);
  }
}

Soft_Cluster_Handler::~Soft_Cluster_Handler() 
{
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
  SP(Cluster) cluster;
  if (clin->size()==1) {
    cluster = (*clin->begin());
    switch (CheckCluster(cluster,true)) {
    case -1:
    case 1:
      msg_Tracking()<<"Potential problem in "<<METHOD<<" : "<<std::endl
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

  Cluster_Iterator cit;
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.), checkaft(0.,0.,0.,0.);
  for (cit=clin->begin();cit!=clin->end();cit++) 
    checkbef += (*cit)->Momentum(); 
#endif

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
  if (count==0) return true;

  double E(-1.);
  bool   tried(false);
  while (!CheckIfAllowed(clin,E)) {
    if (!UpdateTransitions(clin)) {
      if (!TryToEnforceTransition(clin)) {
	if (tried) {
	  if (tried>1) msg_Out()<<"Error in "<<METHOD<<" tried = "<<tried<<"."<<std::endl;
	  return false;
	}
	for (cit=clin->begin();cit!=clin->end();cit++) {
	  cluster = (*cit);
	  if (cluster==NULL || !cluster->Active()) continue;
	  size = CheckCluster(cluster,false,true);
	  switch (size) {
	  case -1:
	    return false;
	  default:
	    count += size;
	  }
	}
	if (count==0) return true;
	tried = true;
      }
      else tried = false;
    }
  }

  if (!ShiftMomenta(clin)) {
    msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		  <<"   Could not shift momenta."<<std::endl
		  <<"   Will possibly lead to retrying the event."<<std::endl;
    return false;
  }

  Particle * part;
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
      part = cluster->GetSelf();
      if (part->Flav()==Flavour(kf_none)) 
	msg_Error()<<"Happens here :"<<part->Flav()<<" vs. "<<(*cluster)[0]<<std::endl;
      blob->AddToOutParticles(part);
#ifdef AHAmomcheck
      checkaft += cluster->Momentum();
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

int Soft_Cluster_Handler::CheckCluster(SP(Cluster) cluster,bool mustdecay,bool lighter)
{
  cluster->clear();

  if (mustdecay) return EnforceDecay(cluster);

  Flavour had1,had2,hadron;
  double decayweight(DecayWeight(cluster,had1,had2));
  double transformweight(TransformWeight(cluster,hadron,lighter));
  if (decayweight==0. && transformweight==0.) return 0;
  if (decayweight<=0. && transformweight!=0.) {
    // no open decay, but open transition (may be forced)
    cluster->push_back(hadron);
    return 1;
  }
  if (decayweight>0. && transformweight<=0.) {
   // no open transition, but open decay
    cluster->push_back(had1);
    cluster->push_back(had2);
    return 2;
  }
  if (decayweight>0. && transformweight>0.) {
    // both transition and decay open
    if (decayweight/(decayweight+transformweight)>ran.Get()) {
      //decay wins
      cluster->push_back(had1);
      cluster->push_back(had2);
      return 2;
    }
    else {
      //transition wins
      if (hadron==Flavour(kf_none)) {
	msg_Error()<<METHOD<<" 2 "<<std::endl;
	abort();
      }
      cluster->push_back(hadron);
      return 1;
    }
  }
  msg_Error()<<"Error in "<<METHOD<<" :"<<std::endl
	     <<(*cluster)<<std::endl
	     <<"   leads to decwt = "<<decayweight
	     <<" and transwt("<<lighter<<") = "<<transformweight<<"."<<std::endl
	     <<"   Trigger new event."<<std::endl;
    
  return -1;
}


int Soft_Cluster_Handler::EnforceDecay(SP(Cluster) cluster) 
{
  Flavour had1,had2;
  double transitweight(DecayWeight(cluster,had1,had2));
  if (transitweight==0.) return 0;
  else if (transitweight>0.) {
    cluster->push_back(had1);
    cluster->push_back(had2);
    return 2;
  }
  else {
    had1          = Flavour(kf_none);
    transitweight = TransformWeight(cluster,had1,true); 
    if (transitweight>0.) {
      cluster->push_back(had1);
      cluster->push_back(Flavour(kf_photon));
      return 2;
    }
  }
  msg_Error()<<"Error in "<<METHOD<<" :"<<std::endl
	     <<"   Did not find any viable decay for "<<std::endl
	     <<(*cluster)<<std::endl
	     <<"   Transition weight = "<<transitweight<<"."<<std::endl
	     <<"   Will trigger new event."<<std::endl;
  return -1;
}

bool Soft_Cluster_Handler::CheckIfAllowed(Cluster_List * clin,double & E) {
  double totmass(0.);
  Vec4D  totmom(0.,0.,0.,0.);
  SP(Cluster) cluster;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
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
  SP(Cluster) cluster;
  Flavour hadron;
  bool    transform(false);
  double  tfwt;
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) {
      hadron = (*cluster)[0];
      tfwt   = TransformWeight(cluster,hadron,true);
      if (tfwt>0.) {
	cluster->clear();
	cluster->push_back(hadron);
	transform = true;
      }
    }
  }
  return transform;
}

bool Soft_Cluster_Handler::TryToEnforceTransition(Cluster_List * clin) {
  bool success(true);
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.),checkaft(0.,0.,0.,0.);
  msg_Out()<<"In "<<METHOD<<" for : "<<std::endl<<(*clin)<<std::endl;
#endif
  size_t size(0);
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  SP(Cluster) cluster;
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
      masses.push_back(cluster->GetTrip()->m_mom.Abs2());
      masses.push_back(cluster->GetAnti()->m_mom.Abs2());
      momenta.push_back(cluster->GetTrip()->m_mom);
      momenta.push_back(cluster->GetAnti()->m_mom);
      size+=2;
    default:
      break;
    }
#ifdef AHAmomcheck
    checkbef += cluster->Momentum();
#endif
  }
  if (!hadpars.AdjustMomenta(size,&momenta.front(),&masses.front())) {
    if (size>1) {
      msg_Tracking()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : "<<std::endl
		    <<"   Could not adjust momenta for : "<<std::endl;
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Debugging()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Debugging()<<" ("<<((**cit)[0])<<" )";
	msg_Debugging()<<" for "<<std::endl<<(**cit)<<std::endl;
      }
      msg_Tracking()<<"   Will possibly lead to retrying the event."<<std::endl;
    }
    return false;
  }
  int pos(0);

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
#ifdef AHAmomcheck
      checkaft += momenta[pos];
#endif
      cluster->GetTrip()->m_mom=momenta[pos++];
      cluster->GetAnti()->m_mom=momenta[pos];
      cluster->Update();		
      break;
    default:
      break;
    }
#ifdef AHAmomcheck
    checkaft += momenta[pos];
#endif
    pos++;
  }
#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields a momentum violation for  "<<size<<" : "<<std::endl
  	     <<"   "<<checkbef<<" - "<<checkaft<<" --> "
  	     <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl
	     <<(*clin)<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum : "
		      <<(checkbef-checkaft).Abs2()<<"("<<size<<") for ."
		      <<std::endl<<(*clin)<<std::endl;
#endif
  return success;  
}


double Soft_Cluster_Handler::DecayWeight(SP(Cluster) cluster,Flavour & had1,Flavour & had2)
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

  double MC(cluster->Mass());
  if (p_doubletransitions->GetHeaviestMass(flpair)<MC+m_offset2) return 0.;
  if (p_doubletransitions->GetLightestMass(flpair)>MC) {
    if (p_singletransitions->MustDesintegrate(cluster,had1,had2)) return 1.;
    if (flpair.first.IsDiQuark() && flpair.second.IsDiQuark())
      msg_Out()<<METHOD<<" for "<<MC<<"["<<flpair.first<<"/"<<flpair.second<<"]"
	       <<" vs. "<<p_doubletransitions->GetLightestMass(flpair)<<std::endl;
    return -1.;
  }
  double totweight(0.),MC2(MC*MC),m1,m2,wt;

  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1  = decit->first.first.HadMass();
    m2  = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * pow((m1+m2)/MC,m_alpha);
      wt *= decit->second;
      totweight += wt;
    }
  }
  if (totweight==0.) {
    msg_Error()<<"Error in "<<METHOD<<" :"<<std::endl
	       <<"   Cluster of mass "<<MC<<" from {"<<flpair.first<<", "
	       <<flpair.second<<"} passed mass conditions,"<<std::endl
	       <<"   but no viable transition found."<<std::endl
	       <<"   Return 0 and hope for the best."<<std::endl;
    throw Return_Value::Retry_Event;
    return 0.;
  }

  had1 = had2 = Flavour(kf_none); 
  double disc(totweight * 0.9999999999*ran.Get());
  for (Double_Transition_Siter decit=dtliter->second->begin();
       decit!=dtliter->second->end();decit++) {
    m1    = decit->first.first.HadMass();
    m2    = decit->first.second.HadMass();
    if (m1+m2<MC) {
      wt  = sqrt((MC2-sqr(m1+m2))*(MC2-sqr(m1-m2))) * pow((m1+m2)/MC,m_alpha);
      wt *= decit->second;
      disc -= wt;
      if (disc<0.) {
	had1 = decit->first.first;
	had2 = decit->first.second;
	break;
      }
    }
  }
  return totweight * 1./(16.*M_PI*MC*MC*MC);
}

double Soft_Cluster_Handler::SelectPT(double pt2max,bool peaked) {
  double pt(0.);
  if (int(m_ptmode)>10 && !peaked) {
    //std::cout<<"   ... isotropic from "<<int(m_ptmode)<<"."<<std::endl;
    double costheta = -1.+2.*ran.Get(), sintheta2 = 1.-sqr(costheta); 
    pt = sqrt(pt2max*sintheta2);
  }
  else {
    //std::cout<<"   ... peaked from "<<int(m_ptmode)<<"."<<std::endl;
    switch (m_ptmode) {
    case (PTdist::sinthet):
    case (PTdist::sinthet_iso):
      pt = sqrt(pt2max)*pow(ran.Get(),m_eta);
      break;
    case (PTdist::alphadist):
    case (PTdist::alphadist_iso):
    default:
      pt = sqrt(p_as->SelectPT(pt2max));
      break;
    }
  }
  return pt;
}

void Soft_Cluster_Handler::FixHHDecay(SP(Cluster) cluster,Blob * blob)
{
#ifdef AHAmomcheck
  Vec4D  checkbef = cluster->Momentum();
#endif
  Flavour had1((*cluster)[0]), had2((*cluster)[1]);

  double M       = cluster->Mass(), M2 = M*M;
  double m12     = sqr(had1.HadMass()), m22 = sqr(had2.HadMass());

  cluster->BoostInCMSAndRotateOnZ();
  double E1      = (M2+m12-m22)/(2.*M);

  bool peaked    = (cluster->GetTrip()->m_info=='L') || (cluster->GetAnti()->m_info=='L');

  //std::cout<<"In "<<METHOD<<" for "<<cluster->GetTrip()->m_info
  //	   <<" "<<cluster->GetAnti()->m_info<<" --> "<<peaked<<"."<<std::endl;

  double pt      = SelectPT(sqr(E1)-m12,peaked);
  double pl1     = sqrt(sqr(E1)-sqr(pt)-m12);
  double cosphi  = cos(2.*M_PI*ran.Get()), sinphi = sqrt(1.-cosphi*cosphi);
  Vec4D  p1      = Vec4D(E1,pt*cosphi,pt*sinphi,pl1);
  Vec4D  p2      = cluster->Momentum()-p1;

  if (m_ana) {
    Histogram* histo((m_histograms.find(std::string("PT_HH")))->second);
    histo->Insert(pt);
    Histogram* histo2((m_histograms.find(std::string("PT2_HH")))->second);
    histo2->Insert(pt*pt);
  }

  Particle * part;
  SP(Cluster) clus;
  clus = new Cluster(p1,had1,false);
  clus->SetPrev(cluster);
  cluster->SetLeft(clus);
#ifdef memchecker
  msg_Debugging()<<"@@@ New cluster "<<clus<<" from "<<METHOD<<"."<<std::endl;
#endif

  clus = new Cluster(p2,had2,false);
  clus->SetPrev(cluster);
  cluster->SetRight(clus);
#ifdef memchecker
  msg_Debugging()<<"@@@ New cluster "<<clus<<" from "<<METHOD<<"."<<std::endl;
#endif

  cluster->RotateAndBoostBack();

  blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
  blob->AddToOutParticles(cluster->GetRight()->GetSelf());


#ifdef AHAmomcheck
  if (dabs((checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2())>1.e-12) {
    msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		  <<"    Four-momentum not conserved: "
		  <<checkbef<<" vs. "<<(cluster->GetLeft()->Momentum()+cluster->GetRight()->Momentum())
		  <<" : "
		  <<(checkbef-cluster->GetLeft()->Momentum()-cluster->GetRight()->Momentum()).Abs2()
		  <<"."<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum:"<<std::endl;
#endif
}


double Soft_Cluster_Handler::TransformWeight(SP(Cluster) cluster,ATOOLS::Flavour & hadron,
					     const bool lighter)
{
  Flavour_Pair fpair;
  fpair.first  = cluster->GetTrip()->m_flav;
  fpair.second = cluster->GetAnti()->m_flav;
  if (fpair.first.IsDiQuark() && fpair.second.IsDiQuark()) return 0.;

  double MC(cluster->Mass());
  if (p_singletransitions->GetHeaviestMass(fpair)+m_offset1<MC) {
    if (!lighter && p_doubletransitions->GetLightestMass(fpair)<MC) return 0.;
  }

  Single_Transition_Miter stiter = p_singletransitions->GetTransitions()->find(fpair);
  if (stiter==p_singletransitions->GetTransitions()->end()) {
    msg_Error()<<"Potential error in  All_Single_Transitions::MustTransit :"<<endl
	       <<"   Did not find any entry for "<<fpair.first<<"/"<<fpair.second
	       <<", mass = "<<cluster->Mass()<<";"<<endl<<(*cluster->GetPrev())
	       <<"   will continue and hope for the best."<<endl;
    return 0.;
  }
  Single_Transition_List * stl(stiter->second);
  Single_Transition_Siter  start(stl->begin()),siter;
  if (lighter && hadron!=Flavour(kf_none)) {
    do {
      if (start->first==hadron) {
	siter = start;
	if ((++siter)!=stl->end()) start++;
	else return -1.;
	break;
      }
      else start++;
    } while (start!=stl->end());
  }

  double minwidth(1.e6),wt,totweight(0.);
  for (Single_Transition_Siter siter=start;siter!=stl->end();siter++) {
    if (!siter->first.IsStable() && (siter->first.Width()<minwidth)) 
      minwidth = siter->first.Width();
  }
  if (minwidth==1.e6) minwidth=1.e-6;
  
  for (siter=start;siter!=stl->end();siter++) {
    if (siter->first!=hadron) {
      if (siter->first.IsStable())
	wt = m_kappa*sqr(MC*minwidth)/
	  //wt = sqr(siter->first.Mass()*minwidth)/
	  (sqr(sqr(MC)-sqr(siter->first.HadMass())) + sqr(siter->first.HadMass()*minwidth));
      else 
	wt = m_kappa*sqr(MC*siter->first.Width())/
	  //wt = sqr(siter->first.Mass()*siter->first.Width())/
	  (sqr(sqr(MC)-sqr(siter->first.HadMass())) + sqr(siter->first.HadMass()*siter->first.Width()));
      wt *= siter->second;
    }
    totweight += wt;
  }

  double disc(totweight * 0.9999999999*ran.Get());
  for (siter=start;siter!=stl->end();siter++) {
    if (siter->first!=hadron) {
      if (siter->first.IsStable())
	wt = m_kappa*sqr(MC*minwidth)/
	  //wt = sqr(siter->first.Mass()*minwidth)/
	  (sqr(sqr(MC)-sqr(siter->first.HadMass())) + sqr(siter->first.HadMass()*minwidth));
      else 
	wt = m_kappa*sqr(MC*siter->first.Width())/
	  //wt = sqr(siter->first.Mass()*siter->first.Width())/
	  (sqr(sqr(MC)-sqr(siter->first.HadMass())) + sqr(siter->first.HadMass()*siter->first.Width()));
      wt *= siter->second;
    }
    disc -= wt;
    if (disc<=0.) {
      hadron = siter->first;
      disc   = wt;
      break;
    }
  }
  return totweight/(2.*MC);;
}

bool Soft_Cluster_Handler::ShiftMomenta(Cluster_List * clin)
{
#ifdef AHAmomcheck
  Vec4D checkbef(0.,0.,0.,0.),checkaft(0.,0.,0.,0.);
#endif
  size_t size(clin->size());
  std::vector<double> masses;
  std::vector<Vec4D>  momenta;
  SP(Cluster) cluster;
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
#ifdef AHAmomcheck
    checkbef += cluster->Momentum();
#endif
  }
  if (!hadpars.AdjustMomenta(size,&momenta.front(),&masses.front())) {
    if (size>1) {
      msg_Tracking()<<"Error in "<<METHOD<<" ("<<size<<" clusters) : "<<std::endl
		    <<"   Could not adjust momenta for : "<<std::endl;
      int i(0);
      for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
	msg_Debugging()<<"Mass/Mom  = "<<masses[i]<<"/"<<momenta[i];
	if ((*cit)->size()==1) msg_Debugging()<<" ("<<((**cit)[0])<<" )";
	msg_Debugging()<<" for "<<std::endl<<(**cit)<<std::endl;
      }
      msg_Tracking()<<"   Will possibly lead to retrying the event."<<std::endl;
    }
    return false;
  }
  int pos(0);
  for (Cluster_Iterator cit=clin->begin();cit!=clin->end();cit++) {
    cluster = (*cit);
    if (cluster->size()==1) cluster->SetFlav((*cluster)[0]);
    else cluster->RescaleMomentum(momenta[pos]);
#ifdef AHAmomcheck
    checkaft += momenta[pos];
#endif
    cluster->SetMomentum(momenta[pos]);
    pos++;
  }
#ifdef AHAmomcheck
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" yields a momentum violation for  "<<size<<" : "<<std::endl
  	     <<"   "<<checkbef<<" - "<<checkaft<<" --> "
  	     <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl;
    msg_Out()<<(*clin)<<std::endl;
  }
  else msg_Debugging()<<METHOD<<" conserves momentum : "
		      <<(checkbef-checkaft).Abs2()<<"("<<size<<")."<<std::endl;
#endif
  return true;
}
