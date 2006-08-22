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
  p_dtransitions(hadpars.GetDoubleTransitions()),
  m_smearparameter(hadpars.Get(string("AngularSmearing")))   
{ }

Return_Value::code Hadron_Part::RedoDecay(Cluster * cluster,Blob * blob,
					  int & mode,Flavour & had1,Flavour & had2)
{
  cout<<"IN "<<METHOD<<" mode = "<<mode<<" for "
      <<cluster->Mass(0)<<" --> "<<had1<<" + "<<had2<<endl;
  p_blob     = blob;
  p_cluster  = cluster;
  m_flavs[0] = had1;
  m_flavs[1] = had2;
  if (m_cht==chtrans::HH_only) return HHDecay();
  switch (mode) {
  case 1: 
  case 2: 
    return CHDecay(mode);
  case 3: return HHDecay(true);
  }
  rvalue.IncError(METHOD);
  return Return_Value::Error;
}

Return_Value::code Hadron_Part::HHDecay(const bool flavs_selected) {
  cout<<"IN "<<METHOD<<endl;
  if (!flavs_selected) {
    if (!p_dtransitions->IsoDecay(p_cluster,m_flavs[0],m_flavs[1]) ||
	m_flavs[0]==Flavour(kf::none) || m_flavs[1]==Flavour(kf::none)) {
      msg.Error()<<"Error in "<<METHOD<<" : "<<endl
		 <<"   Selected an unphysical flavour (none) in decay of"<<endl
		 <<(*p_cluster)<<"   ---> "<<m_flavs[0]<<" + "<<m_flavs[1]<<endl
		 <<"   Will continue and hope for the best."<<endl;
      rvalue.IncError(METHOD);
      return Return_Value::Error;
    }
  }
  if (!CheckHHDecayKinematics()) {  
    msg.Error()<<"Error in "<<METHOD<<" : "<<endl
	       <<"   Selected an unphysical flavour combination in decay of"<<endl
	       <<(*p_cluster)<<"   ---> "<<m_flavs[0]<<" + "<<m_flavs[1]<<endl
	       <<"   Will continue and hope for the best."<<endl;
    rvalue.IncError(METHOD);
    return Return_Value::Error;
  }
  cout<<"   "<<METHOD<<" : "<<p_cluster->Mass(0)<<" --> "
      <<m_flavs[0]<<" + "<<m_flavs[1]<<endl;
  p_cluster->BoostInCMSAndRotateOnZ();
  Vec4D hadmoms[2];
  FillSimpleDecay(hadmoms);
  p_cluster->RotateAndBoostBack(hadmoms[0]);
  p_cluster->RotateAndBoostBack(hadmoms[1]);
  p_cluster->RotateAndBoostBack();
  cout<<"   After back rotate & boost : "<<hadmoms[0]<<" + "<<hadmoms[1]<<endl; 
  FillParticles(hadmoms,3);
  cout<<"   OUT "<<METHOD<<endl<<(*p_blob)<<endl;
  return Return_Value::Success;
}

Return_Value::code Hadron_Part::CHDecay(int & mode) {
  cout<<"IN "<<METHOD<<" mode = "<<mode<<" for "
      <<p_cluster->Mass(0)<<" --> "<<m_flavs[0]<<" + "<<m_flavs[1]<<endl;
  p_cluster->BoostInCMSAndRotateOnZ();
  Vec4D hadmoms[2];
  if (FillCHDecay(hadmoms,mode)) {
    if (mode==1) {
      p_cluster->GetRight()->RescaleMomentum(hadmoms[1]);
      hadmoms[1] = p_cluster->GetRight()->Momentum();
      p_cluster->RotateAndBoostBack(hadmoms[0]);
    }
    else {
      p_cluster->GetLeft()->RescaleMomentum(hadmoms[0]);
      hadmoms[0] = p_cluster->GetLeft()->Momentum();
      p_cluster->RotateAndBoostBack(hadmoms[1]);
    }
    p_cluster->RotateAndBoostBack();
    FillParticles(hadmoms,mode);
    cout<<"   OUT "<<METHOD<<endl<<(*p_blob)<<endl;
    return Return_Value::Success;
  }
  else {
    mode = 3;
    return HHDecay();
  }
}

void Hadron_Part::FillParticles(Vec4D * hadmoms,const int flag) {
  Particle * part;
  for (int i=1;i<3;i++) {
    if (flag&i) {
      part = new Particle(0,m_flavs[i-1],hadmoms[i-1]); 
      part->SetNumber(0);
      part->SetStatus(part_status::active);
      part->SetInfo('P');
      part->SetFinalMass();
      p_blob->AddToOutParticles(part);
    }
  }
}

bool Hadron_Part::CheckHHDecayKinematics()
{
  double m1 = m_flavs[0].Mass(), m2 = m_flavs[1].Mass(), mass = p_cluster->Mass();
  int    next, mode;
  bool   fix(false);
  while (mass<=m1+m2) {
    if (!fix) {
      mode = 0;
      if (ran.Get()>0.5) mode = 1;
    }
    if (mode==0) {
      next = int(p_stransitions->NextLightest(p_cluster->GetLeft(),m_flavs[0]));
      m1   = m_flavs[0].Mass();
    }
    else {
      next = int(p_stransitions->NextLightest(p_cluster->GetRight(),m_flavs[1]));
      m2   = m_flavs[1].Mass();
    }
    switch (next) {
    case (int(Return_Value::Success)) : break;
    case (int(Return_Value::Error))   : 
      if (!fix) {
	mode = 1-mode;
	fix  = true;
	break;
      }
      else return false;  
    }
  }
  return true;
}    


bool Hadron_Part::CheckCHDecayKinematics(Cluster * other)
{
  int mode      = 0;
  Cluster * clu = p_cluster->GetLeft();
  if (clu==other) {
    mode        = 1;
    clu         = p_cluster->GetRight();
  }
  double    m1  = m_flavs[mode].Mass(), m2 = other->Mass(), mass = p_cluster->Mass();
  while (mass<=m1+m2) {
    switch (int(p_stransitions->NextLightest(clu,m_flavs[mode]))) {
    case Return_Value::Success : 
      m1 = m_flavs[mode].Mass(); 
      break;
    case Return_Value::Error : return false;
    }
  }
  return true;
}    

/////////////////////////////////////////////////////////////////////////////////////////////
// Keep_PperpY
/////////////////////////////////////////////////////////////////////////////////////////////
  
Keep_PPerpY::Keep_PPerpY() :
  Hadron_Part()
{
}

void Keep_PPerpY::FillSimpleDecay(ATOOLS::Vec4D* hadmoms) 
{
  double energy   = p_cluster->Momentum()[0];
  double m12      = sqr(m_flavs[0].Mass());
  double m22      = sqr(m_flavs[1].Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);

  double cosphi   = cos(2.*M_PI*ran.Get()), sinphi=sqrt(1.-cosphi*cosphi);
  double pperp    = clu->Momentum().PPerp();
  
  if () {
  }
  else {
    Vec3D  direction(0.,0.,1.);
    double costheta, sintheta, phi;
    if (int(p_cluster->GetLeads())!=0 && m_smearparameter>0.) {
      double norm   = sqr(m_smearparameter);
      do { costheta = 1.-2.*ran.Get(); } while (exp(-(1.-costheta)*norm)<ran.Get());
      sintheta      = sqrt(1.-sqr(costheta));
      direction     = Vec3D(sintheta*cosphi,sintheta*sinphi,costheta);
    }
    hadmoms[0] = Vec4D(energy1,direction*sqrt(sqr(energy1)-m12));
    hadmoms[1] = Vec4D(energy2,(-1.)*direction*sqrt(sqr(energy1)-m12));
  }
}

bool Keep_PPerpY::FillCHDecay(ATOOLS::Vec4D* hadmoms,const int mode)
{
  Cluster * other = (mode==1)?p_cluster->GetRight():p_cluster->GetLeft();
  Cluster * clu   = (mode==1)?p_cluster->GetLeft():p_cluster->GetRight();
  if (!CheckCHDecayKinematics(other)) return false;

  double m2;

  double pperp  = clu->Momentum().PPerp();
  double rap    = clu->Momentum().Y();
  if (mode==1) m2 = sqr(m_flavs[0].Mass());
          else m2 = sqr(m_flavs[1].Mass());
  double mperp  = sqrt(m2+sqr(pperp));
  double E      = mperp*cosh(rap);
  double plong  = mperp*sinh(rap);
  double cosphi = cos(2.*M_PI*ran.Get()), sinphi=sqrt(1.-cosphi*cosphi);

  if (mode==1) {
    hadmoms[0] = Vec4D(E,pperp*cosphi,pperp*sinphi,plong);
    hadmoms[1] = p_cluster->Momentum()-hadmoms[0];
  }
  else {
    hadmoms[1] = Vec4D(E,pperp*cosphi,pperp*sinphi,plong);
    hadmoms[0] = p_cluster->Momentum()-hadmoms[1];
  }
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

void Retain::FillSimpleDecay(ATOOLS::Vec4D * hadmoms) {
  double energy   = p_cluster->Momentum()[0];
  double m12      = sqr(m_flavs[0].Mass());
  double m22      = sqr(m_flavs[1].Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);
  Vec3D  direction(0.,0.,1.);
  double costheta, sintheta, phi;
  if (int(p_cluster->GetLeads())!=0 && m_smearparameter>0.) {
    double norm   = sqr(m_smearparameter);
    do { costheta = 1.-2.*ran.Get(); } while (exp(-(1.-costheta)*norm)<ran.Get());
    sintheta      = sqrt(1.-sqr(costheta));
    phi           = 2.*M_PI*ran.Get();
    direction  = Vec3D(sintheta*sin(phi),sintheta*cos(phi),costheta);
  }
  hadmoms[0] = Vec4D(energy1,direction*sqrt(sqr(energy1)-m12));
  hadmoms[1] = Vec4D(energy2,(-1.)*direction*sqrt(sqr(energy1)-m12));
} 

bool Retain::FillCHDecay(ATOOLS::Vec4D * hadmoms,const int mode) {
  Cluster * other = (mode==1)?p_cluster->GetRight():p_cluster->GetLeft();
  if (!CheckCHDecayKinematics(other)) return false;
  double m12,m22;
  if (mode==1) {
    m12 = sqr(m_flavs[0].Mass());
    m22 = sqr(other->Mass(0));
  }
  else {
    m12 = sqr(other->Mass(0));
    m22 = sqr(m_flavs[1].Mass());
  }

  double energy    = p_cluster->Momentum()[0];
  double energy1 = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2 = (sqr(energy)-m12+m22)/(2.*energy);
  Vec3D  direction(0.,0.,1.);
  hadmoms[0] = Vec4D(energy1,0.,0.,sqrt(sqr(energy1)-m12));
  hadmoms[1] = Vec4D(energy2,0.,0.,-sqrt(sqr(energy1)-m12));
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

void Isotropic::FillSimpleDecay(ATOOLS::Vec4D * hadmoms) {
  double energy   = p_cluster->Momentum()[0];
  double m12      = sqr(m_flavs[0].Mass());
  double m22      = sqr(m_flavs[1].Mass());
  double energy1  = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2  = (sqr(energy)-m12+m22)/(2.*energy);
  Vec3D  direction(0.,0.,1.);
  double costheta, sintheta, phi;
  if (int(p_cluster->GetLeads())==0) {
    costheta        = 1.-2.*ran.Get(); 
    sintheta        = sqrt(1.-sqr(costheta));
    phi             = 2.*M_PI*ran.Get();
  }
  else {
    if (m_smearparameter>0.) {
      double norm   = sqr(m_smearparameter);
      do { costheta = 1.-2.*ran.Get(); } while (exp(-(1.-costheta)*norm)<ran.Get());
      sintheta      = sqrt(1.-sqr(costheta));
      phi           = 2.*M_PI*ran.Get();
    }
    else {
      costheta      = 1.;
      sintheta      = 0.;
      phi           = 0.;
    }
  }
  
  direction  = Vec3D(sintheta*sin(phi),sintheta*cos(phi),costheta);
  hadmoms[0] = Vec4D(energy1,direction*sqrt(sqr(energy1)-m12));
  hadmoms[1] = Vec4D(energy2,(-1.)*direction*sqrt(sqr(energy1)-m12));
} 
   
bool Isotropic::FillCHDecay(ATOOLS::Vec4D * hadmoms,const int mode) {
  Cluster * other = (mode==1)?p_cluster->GetRight():p_cluster->GetLeft();
  if (!CheckCHDecayKinematics(other)) return false;
  double m12,m22;
  if (mode==1) {
    m12 = sqr(m_flavs[0].Mass());
    m22 = sqr(other->Mass(0));
  }
  else {
    m12 = sqr(other->Mass(0));
    m22 = sqr(m_flavs[1].Mass());
  }

  double energy    = p_cluster->Momentum()[0];
  double energy1 = (sqr(energy)+m12-m22)/(2.*energy);
  double energy2 = (sqr(energy)-m12+m22)/(2.*energy);
  Vec3D  direction(0.,0.,1.);
  hadmoms[0] = Vec4D(energy1,0.,0.,sqrt(sqr(energy1)-m12));
  hadmoms[1] = Vec4D(energy2,0.,0.,-sqrt(sqr(energy1)-m12));
} 



