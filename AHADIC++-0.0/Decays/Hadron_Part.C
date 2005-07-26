#include "Hadron_Part.H"
#include "Hadronisation_Parameters.H"
#include "Poincare.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Hadron_Part::Hadron_Part() :
  p_stransitions(hadpars.GetSingleTransitions()),
  p_dtransitions(hadpars.GetDoubleTransitions())   
{ }

void Hadron_Part::RedoDecay(Cluster * cluster,Part_List * pl,int mode,Flavour & had1,Flavour & had2)
{
  if (m_cht==chtrans::HH_only) {
    if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
      TwoHadronDecay(cluster,pl,had1,had2);
      return;
    }
    msg.Error()<<"Error in Hadron_Part::RedoDecay : HH only mode."<<endl
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
    Particle * part(NULL);
    if (mode==2) {
      masses[0] = cluster->GetLeft()->Mass(); 
      masses[1] = had2.Mass();
      switch (CheckDecayKinematics(cluster,cluster->GetLeft(),had2,help)) {
      case 0:
	TwoHadronDecay(cluster,pl,had2,help);
	break;
      case 1:
	hadpars.AdjustMomenta(2,momenta,masses);
	cluster->GetLeft()->RescaleMomentum(momenta[0]); 
	part = new Particle(0,had2,momenta[1]);
	part->SetNumber(0);
	part->SetStatus(1);
	part->SetInfo('P');
	pl->push_back(part);
	break;
      case 2:
	msg.Error()<<"ERROR in Hadron_Part::RedoDecay : "<<endl
		   <<cluster->Mass()<<" -> "<<cluster->GetLeft()->Mass()<<" + "<<had2<<" does not work out .... "<<endl;
	abort();
      default:
	return;
      }
    }
    else {
      masses[0] = had1.Mass();
      masses[1] = cluster->GetRight()->Mass(); 
      switch (CheckDecayKinematics(cluster,cluster->GetRight(),had1,help)) {
      case 0:
	TwoHadronDecay(cluster,pl,had1,help);
	break;
      case 1:
	hadpars.AdjustMomenta(2,momenta,masses);
	cluster->GetRight()->RescaleMomentum(momenta[1]); 
	part = new Particle(0,had1,momenta[0]);
	part->SetNumber(0);
	part->SetStatus(1);
	part->SetInfo('P');
	pl->push_back(part);
	cluster->GetRight()->RescaleMomentum(momenta[1]); 
	break;
      case 2:
	msg.Error()<<"ERROR in Hadron_Part::RedoDecay : "<<endl
		   <<cluster->Mass()<<" -> "<<cluster->GetLeft()->Mass()<<" + "<<had2<<" does not work out .... "<<endl;
	abort();
      default:
	return;
      }
    }
    delete momenta;
    delete masses;
  }
}

bool Hadron_Part::ForcedDecay(Cluster * cluster,Part_List * pl)
{
  Flavour had1, had2;
  if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
    TwoHadronDecay(cluster,pl,had1,had2);
    return true;
  }
  msg.Error()<<"Error in Hadron_Part::ForcedDecay :"<<endl
	     <<"   Could not find a suitable double transition for cluster decay: "
	     <<(*cluster)<<"   will abort the run."<<endl;
  abort();
}

void Hadron_Part::CheckDecayKinematics(Cluster * cluster,Flavour & had1,Flavour & had2)
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


int Hadron_Part::CheckDecayKinematics(Cluster * cluster,Cluster * other,Flavour & had,Flavour & help)
{
  double    m1  = had.Mass(), m2 = other->Mass(), mass = cluster->Mass();
  Cluster * clu = cluster->GetLeft();
  if (clu==other) clu = cluster->GetRight();
  while (mass<=m1+m2) {
    if (!p_stransitions->NextLightest(clu,had)) {
      break;
    }
    else m1 = had.Mass();
  }
  if (mass<=m1+m2) return 1;
 
  help = Flavour(kf::none);
  while (mass<=m1+m2) {
    if (!p_stransitions->NextLightest(other,help)) {
      break;
    }
    else m2 = help.Mass();
  }  
  if (mass>=m1+m2) return 0;
  return 2;
}    

/////////////////////////////////////////////////////////////////////////////////////////////
// Isotropic
/////////////////////////////////////////////////////////////////////////////////////////////

Isotropic::Isotropic() : 
  Hadron_Part(),
  m_smearparameter(hadpars.Get(string("AngularSmearing")))
{
  m_cht    = chtrans::CH_incl; 
  m_hadsel = hadsel::newpair;

  if (m_cht==chtrans::HH_only)   m_hadsel = hadsel::newpair;
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
  Vec3D direction;
  if (int(cluster->GetLeads())>0) {
    if (m_smearparameter>0.) {
      double norm = sqr(m_smearparameter), costheta, sintheta;
      do { costheta = 1.-2.*ran.Get(); } while (exp(-(1.-costheta)*norm)<ran.Get());
      sintheta = sqrt(1.-sqr(costheta));
      double phi      = 2.*M_PI*ran.Get();
      Vec4D help      = Vec4D(1.,sintheta*sin(phi),sintheta*cos(phi),costheta);
      Vec4D dir       = Vec4D(1.,Vec3D(cluster->Momentum(1))/Vec3D(cluster->Momentum(1)).Abs());
      Poincare rot(Vec4D(1.,0.,0.,1.),help);
      rot.Rotate(dir);
      direction       = Vec3D(dir);
    }
    else direction    = Vec3D(cluster->Momentum(1))/Vec3D(cluster->Momentum(1)).Abs();
  }
  else {
    double costheta = 1.-2.*ran.Get(), sintheta = sqrt(1.-sqr(costheta));
    double phi      = 2.*M_PI*ran.Get();
    direction       = Vec3D(sintheta*sin(phi),sintheta*cos(phi),costheta);
  }
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

/////////////////////////////////////////////////////////////////////////////////////////////
// Retain
/////////////////////////////////////////////////////////////////////////////////////////////

Retain::Retain() : 
  Hadron_Part(),
  m_smearparameter(hadpars.Get(string("AngularSmearing")))

{
  m_cht    = chtrans::CH_incl; 
  m_hadsel = hadsel::newpair;

  if (m_cht==chtrans::HH_only)   m_hadsel = hadsel::newpair;
}

void Retain::TwoHadronDecay(Cluster * cluster,Part_List * pl,Flavour & had1,Flavour & had2)
{
  Particle * part;
  cluster->BoostInCMS();
  double energy   = cluster->Momentum()[0];
  double m12      = sqr(had1.Mass());
  double m22      = sqr(had2.Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);
  Vec3D direction;
  if (cluster->GetLeft()!=0) {
    direction = Vec3D(cluster->GetLeft()->Momentum())/(Vec3D(cluster->GetLeft()->Momentum()).Abs());
  }
  else {
    if (int(cluster->GetLeads())>0 && m_smearparameter>0.) {
      double norm = sqr(m_smearparameter), costheta, sintheta;
      do { costheta = 1.-2.*ran.Get(); } while (exp(-(1.-costheta)*norm)<ran.Get());
      sintheta = sqrt(1.-sqr(costheta));
      double phi      = 2.*M_PI*ran.Get();
      Vec4D help      = Vec4D(1.,sintheta*sin(phi),sintheta*cos(phi),costheta);
      Vec4D dir       = Vec4D(1.,Vec3D(cluster->Momentum(1))/Vec3D(cluster->Momentum(1)).Abs());
      Poincare rot(Vec4D(1.,0.,0.,1.),help);
      rot.Rotate(dir);
      direction       = Vec3D(dir);
    }
    else direction    = Vec3D(cluster->Momentum(1))/Vec3D(cluster->Momentum(1)).Abs();
  }
  Vec3D p1        = direction*sqrt(sqr(energy1)-m12);
  Vec3D p2        = (-1.)*p1;
  Vec4D hadmom1   = Vec4D(energy1,p1);
  Vec4D hadmom2   = Vec4D(energy2,p2);
  cluster->BoostBack(hadmom1);
  cluster->BoostBack(hadmom2);
  cluster->BoostBack();
  if (dabs((cluster->Momentum()-hadmom1-hadmom2).Abs2())>1.e-4) {
    msg.Error()<<"Error in Retain::TwoHadronDecay (after boost) : "<<endl
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


