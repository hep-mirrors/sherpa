#include "Cluster_Transformer.H"
#include "Hadronisation_Parameters.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Cluster_Transformer::Cluster_Transformer() :
  m_mode(ctrans::photonemission), 
  m_offset(hadpars.Get(string("Offset"))), m_gammaenergy(0.1), 
  p_stransitions(hadpars.GetSingleTransitions()),
  p_dtransitions(hadpars.GetDoubleTransitions())
{}

Cluster_Transformer::~Cluster_Transformer() {}


Return_Value::code 
Cluster_Transformer::TreatSingleCluster(Cluster *& cluster,Blob * blob) 
{
  Flavour had1, had2;
  bool    decayit(false);
  if (p_dtransitions->MustTransit(cluster,had1,had2)) decayit=true;
  else if (p_stransitions->MustTransit(cluster,had1,m_offset)) {
    double clumass(cluster->Mass()), had1mass(had1.Mass()), had2mass(clumass);
    switch (m_mode) {
      case int(ctrans::photonemission):
	had2 = Flavour(kf::photon);
	while (clumass<had1mass+m_gammaenergy) {
	  switch(int(p_stransitions->NextLightest(cluster,had1))) {
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
	  switch(int(p_stransitions->NextLightest(cluster,had1))) {
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
    cout<<METHOD<<" : decay "<<endl<<(*cluster)<<" --> "<<had1<<" & "<<had2
	<<" ("<<cluster->Mass(0)<<" -> "<<had1.Mass()<<"+"<<had2.Mass()<<")"<<endl;
    DecayCluster(cluster,had1,had2,blob);
//     if (cluster->GetSelf()) { delete cluster->GetSelf(); cluster->SetSelf(NULL); }
//     if (cluster)            { delete cluster; cluster=NULL; }
    return Return_Value::Success;
  }
  return Return_Value::Nothing;
}


Return_Value::code 
Cluster_Transformer::TreatClusterList(Cluster_List * clist,ATOOLS::Blob * blob) {
  Flavour hadron;
  std::map<int,Flavour> hadrons;

  int      number  = clist->size(), i=0;
  Vec4D  * momenta = new Vec4D[number];
  double * masses  = new double[number];
  bool     shiftit = false;
  for (Cluster_Iterator cit=clist->begin();cit!=clist->end();cit++,i++) {
    momenta[i]   = (*cit)->Momentum(0);
    if (p_stransitions->MustTransit((*cit),hadron,m_offset,true)) {
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
      delete masses;
      delete momenta;
      return Return_Value::Error;
    }
    i = 0;
    Particle * part;
    for (Cluster_Iterator cit=clist->begin();cit!=clist->end();i++) {
      if (hadrons.find(i)==hadrons.end()) {
	(*cit)->RescaleMomentum(momenta[i]);  
	blob->AddToOutParticles((*cit)->GetSelf());	
	cit++;
      }
      else {
	part = new Particle(-1,hadrons[i],momenta[i]);
	part->SetNumber();
	part->SetStatus(part_status::active);
	part->SetInfo('P');
	part->SetFinalMass();
	blob->AddToOutParticles(part);
	if ((*cit)->GetSelf()) { delete (*cit)->GetSelf(); (*cit)->SetSelf(NULL); }
	if ((*cit))            { delete (*cit); (*cit)=NULL; }
	cit = clist->erase(cit);
      }
    }
    delete masses;
    delete momenta;
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

  if (dabs((cluster->Momentum()-hadmom1-hadmom2).Abs2())>1.e-4) {
    msg.Error()<<"Error in Isotropic::TwoHadronDecay (after boost) : "<<endl
	       <<"   "<<cluster->Momentum()<<" -> "<<hadmom1<<"+"<<hadmom2<<endl
	       <<"   Sum of moms : "<<cluster->Momentum()-hadmom1-hadmom2<<endl;
  }
  Particle * part = new Particle(-1,had1,hadmom1); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  blob->AddToOutParticles(part);
  part = new Particle(-1,had2,hadmom2); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  blob->AddToOutParticles(part);
}    



