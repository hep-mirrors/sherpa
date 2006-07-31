#include "Cluster_Transformer.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Cluster_Transformer::Cluster_Transformer() :
  m_mode(ctrans::pi0emission), 
  m_offset(hadpars.Get(string("Offset"))), 
  p_transitions(hadpars.GetSingleTransitions())
{}

Cluster_Transformer::~Cluster_Transformer() {}


bool Cluster_Transformer::TreatSingleCluster(Cluster_List * clist,ATOOLS::Blob * blob) {
  Cluster * cluster = *(clist->begin());
  Flavour hadron, extra;
  double  clumass,hadmass,extramass;
  //if (clist->size()==1) cout<<"      Single cluster with "<<cluster->Momentum()<<"."<<endl;  
  if (p_transitions->MustTransit(cluster,hadron,m_offset)) {
    //cout<<"      First cluster with "<<cluster->Momentum()<<"."<<endl;  
    switch (m_mode) {
    case int(ctrans::pi0emission):
      extra     = Flavour(kf::pi);
      clumass   = cluster->Mass();
      hadmass   = hadron.Mass();
      extramass = extra.Mass();
      while (clumass<hadmass+extramass) {
	if (!p_transitions->NextLightest(cluster,hadron)) {
	  //cout<<"Potential error in Cluster_Transformer::TreatSingleCluster."<<std::endl;
	  return true;
	} 
      }
      DecayCluster(cluster,hadron,extra,blob);
      clist->erase(clist->begin());
      return true;
    case int(ctrans::forceddecay):
      msg.Error()<<"Error in  Cluster_Transformer::TreatSingleCluster: "<<endl
		 <<"   Option 'ForcedDecay' not realised yet, will abort."<<endl;
      abort();
    case int(ctrans::nothing):
    default:
      return true;
    }
  }
  return false;
}

bool Cluster_Transformer::TreatSingleCluster(Cluster * cluster,Part_List * plist) {
  Flavour hadron, extra;
  double  clumass,hadmass,extramass;
  //cout<<"      Single cluster with "<<cluster->Momentum()<<" {"
  //  <<cluster->GetFlav(1)<<","<<cluster->GetFlav(2)<<"} "
  //  <<cluster->Momentum().Abs2()<<","<<cluster->Mass()<<"."<<endl;  

  if (p_transitions->MustTransit(cluster,hadron,m_offset)) {
    switch (m_mode) {
    case int(ctrans::pi0emission):
      extra     = Flavour(kf::pi);
      clumass   = cluster->Mass();
      hadmass   = hadron.Mass();
      extramass = extra.Mass();
      while (clumass<hadmass+extramass) {
	if (!p_transitions->NextLightest(cluster,hadron)) {
	  //cout<<"Potential error in Cluster_Transformer::TreatSingleCluster."<<std::endl;
	  return true;
	} 
	hadmass   = hadron.Mass();
      }
      DecayCluster(cluster,hadron,extra,plist);
      return true;
    case int(ctrans::forceddecay):
      msg.Error()<<"Error in  Cluster_Transformer::TreatSingleCluster: "<<endl
		 <<"   Option 'ForcedDecay' not realised yet, will abort."<<endl;
      abort();
    case int(ctrans::nothing):
    default:
      return true;
    }
  }
  return false;
}

void Cluster_Transformer::TreatClusterList(Cluster_List * clist,ATOOLS::Blob * blob) {
  Flavour hadron;
  std::map<int,Flavour> hadrons;

  int      number  = clist->size(), i=0;
  Vec4D  * momenta = new Vec4D[number];
  double * masses  = new double[number];
  bool     shiftit = false;
  //cout<<"Test clusterlist for forced transitions."<<endl;
  for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++,i++) {
    momenta[i]   = (*cit)->Momentum(0);
    //cout<<"Perform check : "<<(*cit)->Mass(0)<<std::endl;
    if (p_transitions->MustTransit((*cit),hadron,m_offset,true)) {
      masses[i]  = hadron.Mass();
      hadrons[i] = hadron;;
      //cout<<"      Cluster in list with "<<(*cit)->Momentum()<<":"<<hadron.Mass()<<endl;  
      shiftit    = true;
    }
    else
      masses[i] = (*cit)->Mass(0);
    //std::cout<<"      Check this : "<<masses[i]<<" "<<(*cit)->Mass(0)<<" "
    //	     <<sqrt(momenta[i].Abs2())<<" "<<sqrt(((*cit)->Momentum(0)).Abs2())<<endl;
  }
  if (shiftit) {
    if (!hadpars.AdjustMomenta(number,momenta,masses)) {
      msg.Error()<<"WARNING in Cluster_Transformer::TreatClusterList :"<<endl
		 <<"   Adjust momenta failed, will continue and hope for the best."<<endl;
    }
    i = 0;
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();i++) {
      if (hadrons.find(i)==hadrons.end()) {
	(*cit)->RescaleMomentum(momenta[i]);  
	cit++;
      }
      else {
	Particle * part = new Particle(0,hadrons[i],momenta[i]);
	part->SetNumber(0);
	part->SetStatus(part_status::active);
	part->SetInfo('P');
	part->SetFinalMass();
	blob->AddToOutParticles(part);
	cit = clist->erase(cit);
      }
    }
  }
  delete masses;
  delete momenta;
}

void Cluster_Transformer::DecayCluster(Cluster * cluster,Flavour & had1,Flavour & had2,
				       Blob * blob)
{
  //cout<<"         Check C->HH (1): "<<cluster->Momentum()<<" -> "<<had1<<"/"<<had2<<endl;
  cluster->BoostInCMS();
  double energy   = cluster->Momentum()[0];
  double m12      = sqr(had1.Mass());
  double m22      = sqr(had2.Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);
  double costheta = 1.-2.*ran.Get(), sintheta = sqrt(1.-sqr(costheta));
  double phi      = 2.*M_PI*ran.Get();
  Vec3D direction = Vec3D(sintheta*sin(phi),sintheta*cos(phi),costheta);
  Vec3D p1        = direction*sqrt(sqr(energy1)-m12);
  Vec3D p2        = (-1.)*p1;
  Vec4D hadmom1   = Vec4D(energy1,p1);
  Vec4D hadmom2   = Vec4D(energy2,p2);
  cluster->BoostBack(hadmom1);
  cluster->BoostBack(hadmom2);
  cluster->BoostBack();

  //cout<<"         Check C->HH (2): "<<cluster->Momentum()
  //<<endl<<"   -> "<<hadmom1<<" "<<hadmom2<<endl;
  
  if (dabs((cluster->Momentum()-hadmom1-hadmom2).Abs2())>1.e-4) {
    msg.Error()<<"Error in Isotropic::TwoHadronDecay (after boost) : "<<endl
	       <<"   "<<cluster->Momentum()<<" -> "<<hadmom1<<"+"<<hadmom2<<endl
	       <<"   Sum of moms : "<<cluster->Momentum()-hadmom1-hadmom2<<endl;
  }
  Particle * part = new Particle(0,had1,hadmom1); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  blob->AddToOutParticles(part);
  part = new Particle(0,had2,hadmom2); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  blob->AddToOutParticles(part);

  //cout<<"         Filling blob from Cluster_Transformer : "<<(*blob)<<endl;
}    


void Cluster_Transformer::DecayCluster(Cluster * cluster,Flavour & had1,Flavour & had2,Part_List * plist)
{
  //cout<<"         Check C->HH (1): "<<cluster->Momentum()<<endl;
  cluster->BoostInCMS();
  double energy   = cluster->Momentum()[0];
  double m12      = sqr(had1.Mass());
  double m22      = sqr(had2.Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);
  double costheta = 1.-2.*ran.Get(), sintheta = sqrt(1.-sqr(costheta));
  double phi      = 2.*M_PI*ran.Get();
  Vec3D direction = Vec3D(sintheta*sin(phi),sintheta*cos(phi),costheta);
  Vec3D p1        = direction*sqrt(sqr(energy1)-m12);
  Vec3D p2        = (-1.)*p1;
  Vec4D hadmom1   = Vec4D(energy1,p1);
  Vec4D hadmom2   = Vec4D(energy2,p2);
  cluster->BoostBack(hadmom1);
  cluster->BoostBack(hadmom2);
  cluster->BoostBack();

  //cout<<"         Check C->HH (2): "<<cluster->Momentum()<<endl<<"   -> "<<hadmom1<<" "<<hadmom2<<endl;
  
  if (dabs((cluster->Momentum()-hadmom1-hadmom2).Abs2())>1.e-4) {
    msg.Error()<<"Error in Isotropic::TwoHadronDecay (after boost) : "<<endl
	       <<"   "<<cluster->Momentum()<<" -> "<<hadmom1<<"+"<<hadmom2<<endl
	       <<"   Sum of moms : "<<cluster->Momentum()-hadmom1-hadmom2<<endl;
  }
  Particle * part = new Particle(0,had1,hadmom1); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  plist->push_back(part);
  part = new Particle(0,had2,hadmom2); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  plist->push_back(part);

  //cout<<"         Filling blob from Cluster_Transformer : "<<(*blob)<<endl;
}    

