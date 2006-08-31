#include "Cluster_Transformer.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Cluster_Transformer::Cluster_Transformer() :
  m_mode(ctrans::photonemission), 
  m_offset(hadpars.Get(string("Offset_C->H"))), m_gammaenergy(0.1), 
  p_transitions(hadpars.GetSingleTransitions()),
  p_popper(hadpars.GetPopper())
{}

Cluster_Transformer::~Cluster_Transformer() {}


Return_Value::code 
Cluster_Transformer::TreatSingleCluster(Cluster *& cluster,Blob * blob) 
{
  Flavour had1, had2;
  bool    decayit(false);
  if (p_transitions->MustDesintegrate(cluster,had1,had2)) {
    decayit = true;
  }
  else if (p_transitions->MustTransit(cluster,had1,m_offset)) {
    double clumass(cluster->Mass()), had1mass(had1.Mass()), had2mass(clumass);
    switch (m_mode) {
    case int(ctrans::photonemission):
      had2 = Flavour(kf::photon);
      while (clumass<had1mass+m_gammaenergy) {
	switch(int(p_transitions->NextLightest(cluster,had1))) {
	case (int(Return_Value::Success)) :
	  break;
	case (int(Return_Value::Error))   :
	case (int(Return_Value::Warning)) :
	  return Return_Value::Warning;
	}
	had1mass   = had1.Mass();
      }
      decayit = true;
      break;
    case int(ctrans::pi0emission):
      had2     = Flavour(kf::pi);
      had2mass = had2.Mass();
      while (clumass<had1mass+had2mass) {
	switch(int(p_transitions->NextLightest(cluster,had1))) {
	case (int(Return_Value::Success)) :
	  break;
	case (int(Return_Value::Error))   :
	case (int(Return_Value::Warning)) :
	  return Return_Value::Warning;
	}
	had1mass   = had1.Mass();
      }
      decayit = true;
      break;
    default:
      return Return_Value::Warning;
    }
  }
  
  if (decayit) {
    DecayCluster(cluster,had1,had2,blob);
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}


Return_Value::code 
Cluster_Transformer::TreatClusterList(Cluster_List * clist,ATOOLS::Blob * blob) {

  Cluster_Iterator cit=clist->begin();
  Flavour hadron, hadron2;
  while (cit!=clist->end()) {
    if (p_transitions->MustDesintegrate((*cit),hadron,hadron2)) {
      DecayCluster((*cit),hadron,hadron2,blob);
      if ((*cit)->GetSelf()) { 	  
	control::s_AHAparticles--;
	delete (*cit)->GetSelf(); (*cit)->SetSelf(NULL); 
      }
      if ((*cit)) { 
	delete (*cit); (*cit)=NULL; 
      }
      cit = clist->erase(cit);
    }
    else cit++;
  }

  std::map<int,Flavour> hadrons;
  int    number  = clist->size(), i=0;
  Vec4D  momenta[number];
  double masses[number];
  bool   shiftit = false;

  for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++,i++) {
    momenta[i]   = (*cit)->Momentum(0);
    if (p_transitions->MustTransit((*cit),hadron,m_offset,true)) {
      masses[i]  = hadron.Mass();
      hadrons[i] = hadron;;
      shiftit    = true;
    }
    else masses[i] = (*cit)->Mass(0);
  }
  if (shiftit) {
    if (!hadpars.AdjustMomenta(number,momenta,masses)) {
      msg.Error()<<"WARNING in Cluster_Transformer::TreatClusterList :"<<endl
		 <<"   Adjust momenta failed, will continue and hope for the best."<<endl;
      return Return_Value::Error;
    }

    i = 0;
    Particle * part;
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();i++) {
      if (hadrons.find(i)==hadrons.end()) {
	(*cit)->RescaleMomentum(momenta[i]);  
	(*cit)->GetSelf()->SetMomentum((*cit)->Momentum());
	blob->AddToOutParticles((*cit)->GetSelf());	
	cit++;
      }
      else {
	part = new Particle(-1,hadrons[i],momenta[i]);
	part->SetNumber();
	part->SetStatus(part_status::active);
	part->SetInfo('P');
	part->SetFinalMass();
	control::s_AHAparticles++;
	blob->AddToOutParticles(part);
	if ((*cit)->GetSelf()) { 	  
	  control::s_AHAparticles--;
	  delete (*cit)->GetSelf(); (*cit)->SetSelf(NULL); 
	}
	if ((*cit)) { 
	  delete (*cit); (*cit)=NULL; 
	}
	cit = clist->erase(cit);
      }
    }
  }
  else {
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++) 
      blob->AddToOutParticles((*cit)->GetSelf());	
  }
  return Return_Value::Success;
}

void Cluster_Transformer::DecayCluster(Cluster * cluster,Flavour & had1,Flavour & had2,
				       Blob * blob)
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

  cluster->RotateAndBoostBack(p1);
  cluster->RotateAndBoostBack(p2);

  Particle * part;
  part = new Particle(-1,had1,p1);
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had1.PSMass());
  blob->AddToOutParticles(part);
  control::s_AHAparticles++;
  part = new Particle(-1,had2,p2);
  part->SetNumber();
  part->SetInfo('P');
  part->SetStatus(part_status::active);
  part->SetFinalMass(had2.PSMass());
  blob->AddToOutParticles(part);
  control::s_AHAparticles++;
}    



