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

Return_Value::code  Hadron_Part::RedoDecay(Cluster * cluster,Part_List * pl,
					   int mode,Flavour & had1,Flavour & had2)
{
  if (m_cht==chtrans::HH_only) return HHDecay(cluster,pl,had1,had2);
  if (mode==3)                 return CHDecayBoth(cluster,pl,had1,had2);

  if ((mode==1 && had1==Flavour(kf::none)) || 
      (mode==2 && had2==Flavour(kf::none))) {
    msg.Error()<<"Error in "<<METHOD<<"(mode = "<<mode<<") : "<<endl
	       <<"   Selected an unphysical flavour (none) in decay of"<<endl
	       <<(*cluster)<<"   ---> "<<had1<<" + "<<had2<<endl
	       <<"   Will continue and hope for the best."<<endl;
    rvalue.IncError(METHOD);
    return Return_Value::Error;
  }
  else if (mode==1)            return CHDecayOne(cluster,pl,cluster->GetRight(),had1);
  else if (mode==2)            return CHDecayOne(cluster,pl,cluster->GetLeft(),had2);
  rvalue.IncError(METHOD);
  //cout<<METHOD<<" 1: Error."<<endl;
  return Return_Value::Error;
}

Return_Value::code Hadron_Part::ForcedDecay(Cluster * cluster,Part_List * pl)
{
  Flavour had1, had2;
  return HHDecay(cluster,pl,had1,had2);
}

Return_Value::code Hadron_Part::CheckDecayKinematics(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  Flavour   had;
  Cluster * clu;
  bool      mode;
  double    m1 = had1.Mass(), m2 = had2.Mass(), mass = cluster->Mass();
  while (mass<=m1+m2) {
    if (ran.Get()>0.5) { had = had2; clu=cluster->GetRight(); mode = 1; }
                  else { had = had1; clu=cluster->GetLeft();  mode = 0; }
    switch (int(p_stransitions->NextLightest(clu,had))) {
    case (int(Return_Value::Success)) :  
      if (mode) { had2 = had; m2 = had2.Mass(); } 
           else { had1 = had; m1 = had1.Mass(); } 
      break;
    case (int(Return_Value::Error)) :
      return Return_Value::Error;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
      break;
    }
  }
  return Return_Value::Success;
}    


Return_Value::code  Hadron_Part::CheckDecayKinematics(Cluster * cluster,Cluster * other,
						      Flavour & had,Flavour & help)
{
  //cout<<METHOD<<endl;
  double    m1  = had.Mass(), m2 = other->Mass(), mass = cluster->Mass();
  Cluster * clu(other==cluster->GetLeft()?cluster->GetRight():cluster->GetLeft());
  bool      exit(false);
  while (mass<=m1+m2 && !exit) {
    switch (int(p_stransitions->NextLightest(clu,had))) {
    case Return_Value::Success : 
      m1 = had.Mass(); break;
    case Return_Value::Error : return Return_Value::Error; 
    case Return_Value::Warning : 
      exit=true; break;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
      break;
    }
  }
  if (mass>=m1+m2) return Return_Value::Success; 
  exit = false;
  while (mass<=m1+m2 && !exit) {
    switch (int(p_stransitions->NextLightest(other,help))) {
    case Return_Value::Success : 
      m2 = help.Mass(); break;
    case Return_Value::Error : return Return_Value::Error; 
    case Return_Value::Warning : 
      exit=true; break;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
      break;
    }
  }  
  if (mass>=m1+m2) return Return_Value::Nothing; 
  return Return_Value::Error; 
}    


Return_Value::code Hadron_Part::HHDecay(Cluster * cluster,Part_List * pl,Flavour & had1,Flavour & had2) {
  //cout<<METHOD<<endl;
  if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
    if (had1==Flavour(kf::none) || had2==Flavour(kf::none)) {
      msg.Error()<<"Error in "<<METHOD<<" : "<<endl
		 <<"   Selected an unphysical flavour (none) in decay of"<<endl
		 <<(*cluster)<<"   ---> "<<had1<<" + "<<had2<<endl
		 <<"   Will continue and hope for the best."<<endl;
      rvalue.IncError(METHOD);
      //cout<<METHOD<<" 1: Error."<<endl;
      return Return_Value::Error;
    }
    TwoHadronDecay(cluster,pl,had1,had2);
    //cout<<METHOD<<" 1: Success."<<endl;
    return Return_Value::Success;
  }
  msg.Error()<<"Error in "<<METHOD<<" : "<<endl
	     <<"   Could not find a suitable double transition for cluster decay: "
	     <<(*cluster)<<"   will continue and hope for the best."<<endl;
  //cout<<METHOD<<" 2: Error."<<endl;
  rvalue.IncError(METHOD);
  return Return_Value::Error;
}

Return_Value::code Hadron_Part::CHDecayOne(Cluster * cluster,Part_List * pl,
					   Cluster * remain,Flavour & had) {
  //cout<<METHOD<<endl;
  Vec4D  * momenta = new Vec4D[2];
  double * masses  = new double[2];
  Flavour help = Flavour(kf::none);
  Particle * part(NULL);
  Cluster * partner(cluster->GetLeft()==remain?cluster->GetRight():cluster->GetLeft());

  momenta[0]   = remain->Momentum();
  momenta[1]   = partner->Momentum();
  masses[0]    = remain->Mass(); 
  masses[1]    = had.PSMass();
  switch (int(CheckDecayKinematics(cluster,remain,had,help))) {
  case int(Return_Value::Success):
    hadpars.AdjustMomenta(2,momenta,masses);
    remain->RescaleMomentum(momenta[0]); 
    part = new Particle(0,had,momenta[1]);
    part->SetNumber(0);
    part->SetStatus(part_status::active);
    part->SetInfo('P');
    part->SetFinalMass();
    pl->push_back(part);
    break;
  case int(Return_Value::Nothing):
    TwoHadronDecay(cluster,pl,had,help);
    break;
  case int(Return_Value::Error):
    delete momenta;
    delete masses;
    msg.Error()<<"ERROR in Hadron_Part::RedoDecay : "<<endl
	       <<cluster->Mass()<<" -> "
	       <<cluster->GetLeft()->Mass()<<" + "<<had
	       <<" does not work out .... "<<endl;
    rvalue.IncError(METHOD);
    //cout<<METHOD<<" 1: Error."<<endl;
    return Return_Value::Error;
  default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
    break;
  }
  delete momenta;
  delete masses;
  //cout<<METHOD<<" 1: Success."<<endl;
  return Return_Value::Success;
}


Return_Value::code Hadron_Part::CHDecayBoth(Cluster * cluster,Part_List * pl,
					    Flavour & had1,Flavour & had2) {
  //cout<<METHOD<<endl;
  if (m_hadsel==hadsel::newpair) {
    if (p_dtransitions->IsoDecay(cluster,had1,had2)) { 
      if (had1==Flavour(kf::none) || had2==Flavour(kf::none)) {
	msg.Error()<<"Error in "<<METHOD<<" : "<<endl
		   <<"   Selected an unphysical flavour (none) in decay of"<<endl
		   <<(*cluster)<<"   ---> "<<had1<<" + "<<had2<<endl
		   <<"   Will continue and hope for the best."<<endl;
	rvalue.IncError(METHOD);
	//cout<<METHOD<<" 1: Error."<<endl;
	return Return_Value::Error;
      }
    }
  }
  else {
    switch (int(CheckDecayKinematics(cluster,had1,had2))) {
    case int(Return_Value::Success) :
    case int(Return_Value::Nothing) :
      break;
    case int(Return_Value::Error) :
      rvalue.IncError(METHOD);
      //cout<<METHOD<<" 2: Error."<<endl;
      return Return_Value::Error;
    default:
      msg.Error()<<"Error in "<<METHOD<<": "<<endl
		 <<"   Unknown return value."<<endl;
      abort();
    }
  }
  TwoHadronDecay(cluster,pl,had1,had2);
  //cout<<METHOD<<" 1: Success."<<endl;
  return Return_Value::Success;
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
    rvalue.IncWarning(METHOD);
  }
  part = new Particle(0,had1,hadmom1); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  pl->push_back(part);
  //cout<<"TwoHadronDecay : Add "<<had1<<" to pl"<<endl;
  part = new Particle(0,had2,hadmom2); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  pl->push_back(part);
  //cout<<"TwoHadronDecay : Add "<<had2<<" to pl"<<endl;
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
    rvalue.IncWarning(METHOD);
  }
  part = new Particle(0,had1,hadmom1); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  //cout<<"TwoHadronDecay : Add "<<had1<<" to pl"<<endl;
  pl->push_back(part);
  part = new Particle(0,had2,hadmom2); 
  part->SetNumber(0);
  part->SetStatus(part_status::active);
  part->SetInfo('P');
  part->SetFinalMass();
  pl->push_back(part);
  //cout<<"TwoHadronDecay : Add "<<had2<<" to pl"<<endl;
}    


