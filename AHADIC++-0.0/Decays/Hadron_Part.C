#include "Hadron_Part.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Hadron_Part::Hadron_Part() :
  p_stransitions(hadpars.GetSingleTransitions()),
  p_dtransitions(hadpars.GetDoubleTransitions())   
{ }

Isotropic::Isotropic() : Hadron_Part()
{
  m_cht    = chtrans::CH_incl; 
  m_hadsel = hadsel::newpair;

  if (m_cht==chtrans::HH_only)   m_hadsel = hadsel::newpair;
}

void Isotropic::RedoDecay(Cluster * cluster,Part_List * pl,int mode,Flavour & had1,Flavour & had2)
{
  if (m_cht==chtrans::HH_only) {
    if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
      TwoHadronDecay(cluster,pl,had1,had2);
      return;
    }
    msg.Error()<<"Error in Isotropic::RedoDecay : HH only mode."<<endl
	       <<"   Could not find a suitable double transition for cluster decay: "
	       <<(*cluster)<<"   will abort the run."<<endl;
    abort();
  }
  if (mode==3) {
    if (m_hadsel==hadsel::newpair) {
      if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
      }
    }
    else CheckDecayKinematics(cluster,had1,had2);

    TwoHadronDecay(cluster,pl,had1,had2);
  }
  else {
    Vec4D  * momenta = new Vec4D[2];
    double * masses  = new double[2];
    momenta[0]  = cluster->GetLeft()->Momentum();
    momenta[1]  = cluster->GetRight()->Momentum();
    Flavour help;
    if (mode==2) {
      masses[0] = cluster->GetLeft()->Mass(); 
      masses[1] = had2.Mass();
      if (CheckDecayKinematics(cluster,cluster->GetLeft(),had2,help)) {
	TwoHadronDecay(cluster,pl,had2,help);
	return;
      }
    }
    else {
      masses[0] = had1.Mass();
      masses[1] = cluster->GetRight()->Mass(); 
      if (CheckDecayKinematics(cluster,cluster->GetRight(),had1,help)) {
	TwoHadronDecay(cluster,pl,had1,had2);
	return;
      }
    }
    hadpars.AdjustMomenta(2,momenta,masses);
    if (mode==2) {
      cluster->GetLeft()->RescaleMomentum(momenta[0]); 
      Particle * part = new Particle(0,had2,momenta[1]);
      part->SetNumber(0);
      part->SetStatus(1);
      part->SetInfo('P');
      pl->push_back(part);
    }
    else {
      Particle * part = new Particle(0,had1,momenta[0]);
      part->SetNumber(0);
      part->SetStatus(1);
      part->SetInfo('P');
      pl->push_back(part);
      cluster->GetRight()->RescaleMomentum(momenta[1]); 
    }
    delete momenta;
    delete masses;
  }
}

bool Isotropic::ForcedDecay(Cluster * cluster,Part_List * pl)
{
  Flavour had1, had2;
  if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
    TwoHadronDecay(cluster,pl,had1,had2);
    return true;
  }
  msg.Error()<<"Error in Isotropic::ForcedDecay :"<<endl
	     <<"   Could not find a suitable double transition for cluster decay: "
	     <<(*cluster)<<"   will abort the run."<<endl;
  abort();
}

void Isotropic::TwoHadronDecay(Cluster * cluster,Part_List * pl,Flavour & had1,Flavour & had2)
{
  Particle * part;
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
  part = new Particle(0,had1,hadmom1); 
  part->SetNumber(0);
  part->SetStatus(1);
  part->SetInfo('P');
  pl->push_back(part);
  part = new Particle(0,had2,hadmom2); 
  part->SetNumber(0);
  part->SetStatus(1);
  part->SetInfo('P');
  pl->push_back(part);
}    



void Isotropic::CheckDecayKinematics(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  Flavour   had;
  Cluster * clu;
  bool      mode;
  double    m1 = had1.Mass(), m2 = had2.Mass(), mass = cluster->Mass();
  while (mass<=m1+m2) {
    if (ran.Get()>0.5) { had = had2; clu=cluster->GetRight(); mode = 1; }
                  else { had = had1; clu=cluster->GetLeft();  mode = 0; }
    if (p_stransitions->NextLightest(clu,had)) {
      if (mode) { had2 = had; m2 = had2.Mass(); } 
           else { had1 = had; m1 = had1.Mass(); } 
    }
  }
  return;
}    


bool Isotropic::CheckDecayKinematics(Cluster * cluster,Cluster * other,Flavour & had,Flavour & help)
{
  double    m1  = had.Mass(), m2 = other->Mass(), mass = cluster->Mass();
  Cluster * clu = cluster->GetLeft();
  if (clu==other) clu = cluster->GetRight();
  while (mass<=m1+m2) {
    if (!p_stransitions->NextLightest(clu,had)) {
      break;
    }
    m1 = had.Mass();
  }
  if (mass>=m1+m2) false;

  help = Flavour(kf::none);
  while (mass<=m1+m2) {
    if (!p_stransitions->NextLightest(other,help)) {
    }
    return true;
  }  
}    

