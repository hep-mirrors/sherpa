#include "Cluster_Decay_Mode.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;

Cluster_Decay_Mode::Cluster_Decay_Mode() :
  p_transformer(NULL), p_hadrons(NULL)
{}

Cluster_Decay_Mode::Cluster_Decay_Mode(Clusters_2_Hadrons * _transformer) :
  p_transformer(_transformer), p_hadrons(new std::deque<ATOOLS::Particle *>)
{
  m_totalwt    = 0.;
  for (FlavCCMap_Iterator fliter=hadpars.GetConstituents()->CCMap.begin();
       fliter!=hadpars.GetConstituents()->CCMap.end();fliter++) {
    if ((fliter->first.IsQuark() || fliter->first.IsDiQuark()) && 
	fliter->second->Mass()<hadpars.Get("Max_Prod_Mass")) { 
      m_dicelist.insert(std::make_pair(fliter->first,fliter->second));
      m_totalwt+=fliter->second->TotWeight();
    }
  }
}

Cluster_Decay_Mode::~Cluster_Decay_Mode() 
{
  Reset();
  if (p_hadrons)  { delete p_hadrons;  p_hadrons  = NULL; } 
  if (!m_dicelist.empty()) {
    for (FlavCCMap_Iterator fliter=m_dicelist.begin();fliter!=m_dicelist.end();fliter++) {  
      if (fliter->second) {delete fliter->second; fliter->second=NULL; }
    }
    m_dicelist.clear();
  }
}



Simple_Cluster_Fission::Simple_Cluster_Fission(Clusters_2_Hadrons * _transformer) :
  Cluster_Decay_Mode(_transformer), m_Q(1.), p_masses(new double[2])
{ 
}

Simple_Cluster_Fission::~Simple_Cluster_Fission() 
{
  if (p_masses) { delete [] p_masses; p_masses=NULL; }
}

void Simple_Cluster_Fission::Decay(Cluster * cluster) 
{
  cluster->BoostInCMS();
  ConstructTestClusters(cluster);
  cluster->BoostBack();
  switch (p_transformer->Transition(cluster)) {
  case 3: // Two hadrons
    IsotropicDecay(cluster); 
    break;
  case 2: // One hadron (B)   
    RearrangeMomenta(cluster);
    Decay(cluster->p_left);
    break; 
  case 1: // One hadron (A)
    RearrangeMomenta(cluster);
    Decay(cluster->p_right);
    break;
  case 0: // No hadrons  
    RearrangeMomenta(cluster);
    Decay(cluster->p_left);
    Decay(cluster->p_right);
    break;
  }
}

void Simple_Cluster_Fission::ConstructTestClusters(Cluster * cluster) 
{
  double mass      = sqrt(cluster->Mass2());
  double mass1     = hadpars.GetConstituents()->Mass(cluster->FlavPair()->first);
  double mass2     = hadpars.GetConstituents()->Mass(cluster->FlavPair()->second);
  SelectFlavour(mass-mass1-mass2);
    

  for (int i=0;i<2;i++) m_inmoms[i] = cluster->Momentum(i);
  m_outmoms[0]     = (1.-m_Q/mass)*m_inmoms[0];
  m_outmoms[1]     = m_Q/mass*m_inmoms[1];
  m_outmoms[2]     = m_Q/mass*m_inmoms[0];
  m_outmoms[3]     = (1.-m_Q/mass)*m_inmoms[1];

  p_masses[0]      = sqrt((m_outmoms[0]+m_outmoms[1]).Abs2());
  p_masses[1]      = sqrt((m_outmoms[2]+m_outmoms[3]).Abs2());

  Particle * part1 = NULL, * part2 = NULL;
  part1            = new Particle(cluster->GetParticle(0)->Number(), 
				  cluster->GetParticle(0)->Flav(),m_outmoms[0], 
				  cluster->GetParticle(0)->Info());
  part2            = new Particle(-1,m_flav.Bar(),m_outmoms[1],'c');
  part2->SetNumber(0);
  cluster->p_left  = new Cluster(part1,part2);
  part1            = new Particle(-1,m_flav,m_outmoms[2],'c');
  part1->SetNumber(0);
  part2            = new Particle(cluster->GetParticle(1)->Number(), 
				  cluster->GetParticle(1)->Flav(),m_outmoms[3], 
				  cluster->GetParticle(1)->Info());
  cluster->p_right = new Cluster(part1,part2);

  Vec4D testmom = cluster->Momentum()-cluster->p_left->Momentum()-cluster->p_right->Momentum();
}

void Simple_Cluster_Fission::SelectFlavour(double maxmass) 
{
  double wt;
  FlavCCMap_Iterator fliter;
  for (;;) {
    wt = ran.Get()*m_totalwt;
    for (fliter=m_dicelist.begin();fliter!=m_dicelist.end();fliter++) {
      if (2.001*fliter->second->Mass()<maxmass) wt -= fliter->second->TotWeight();
      if (wt<0.) { 
	m_flav = fliter->first; 
	if (m_flav.IsDiQuark()) m_flav = m_flav.Bar();
	return; 
      }
    }
  }
}

void Simple_Cluster_Fission::RearrangeMomenta(Cluster * cluster)
{
  double M2, mm12, mm22, pabs;
  Vec4D p1, p2;
  Vec3D dir;
  Cluster * cl;
  Vec4D testmom;
  if ((!cluster->p_left->FlavPair()->first.IsHadron()) &&
      (!cluster->p_right->FlavPair()->first.IsHadron())) {
    for (int i=0;i<2;i++) {
      if (i==0) cl = cluster->p_left;
      if (i==1) cl = cluster->p_right;
      M2   = cl->Mass2(); 
      mm12 = sqr(hadpars.GetConstituents()->Mass(cl->FlavPair()->first)); 
      mm22 = sqr(hadpars.GetConstituents()->Mass(cl->FlavPair()->second)); 
      if (sqrt(M2)<sqrt(mm12+mm22)) {
	msg.Error()<<"Error in Simple_Cluster_Fission::RearrangeMomenta."<<std::endl
		   <<"   Funny masses : "<<M2<<"  <  "<<mm12<<" + "<<mm22<<std::endl;
	abort();
      }
      cl->BoostInCMS();
      dir  = Vec3D(cl->Momentum(0)); dir = dir/dir.Abs();
      pabs = sqrt(sqr(M2-mm12-mm22)-4.*mm12*mm22)/(2.*sqrt(M2));
      p1   = Vec4D((M2-mm22+mm12)/(2.*sqrt(M2)),pabs*dir);
      p2   = Vec4D((M2+mm22-mm12)/(2.*sqrt(M2)),-pabs*dir);
      cl->GetParticle(0)->SetMomentum(p1);
      cl->GetParticle(1)->SetMomentum(p2);
      cl->BoostBack();
    }
    testmom = cluster->Momentum()-cluster->p_left->Momentum()-cluster->p_right->Momentum();
  }
  else if ((cluster->p_left->FlavPair()->first.IsHadron()) &&
	   (!cluster->p_right->FlavPair()->first.IsHadron())) {
    Vec4D help = cluster->Momentum(), help1 = cluster->p_left->Momentum(), help2 = cluster->p_right->Momentum();

    cluster->BoostInCMS();
    M2   = cluster->Mass2();
    mm12 = sqr(cluster->p_left->FlavPair()->first.Mass());
    mm22 = cluster->p_right->Mass2();
    dir  = Vec3D(cluster->p_left->Momentum()); dir = dir/dir.Abs();
    pabs = sqrt(sqr(M2-mm12-mm22)-4.*mm12*mm22)/(2.*sqrt(M2));
    p1   = Vec4D((M2-mm22+mm12)/(2.*sqrt(M2)),pabs*dir);
    p2   = Vec4D((M2+mm22-mm12)/(2.*sqrt(M2)),-pabs*dir);
    cluster->p_right->RescaleMomentum(p2);
    cluster->BoostBack(p1); cluster->BoostBack();
    Particle * h1 = new Particle(-1,cluster->p_left->FlavPair()->first,p1,'P');
    h1->SetNumber(0);
    p_hadrons->push_back(h1);
    testmom = cluster->Momentum()-p1-cluster->p_right->Momentum();
  }
  else if (!(cluster->p_left->FlavPair()->first.IsHadron()) &&
	   (cluster->p_right->FlavPair()->first.IsHadron())) {
    Vec4D help = cluster->Momentum(), help1 = cluster->p_left->Momentum(), help2 = cluster->p_right->Momentum();
    cluster->BoostInCMS();
    M2   = cluster->Mass2();
    mm12 = cluster->p_left->Mass2();
    mm22 = sqr(cluster->p_right->FlavPair()->first.Mass());
    dir  = Vec3D(cluster->p_left->Momentum()); dir = dir/dir.Abs();
    pabs = sqrt(sqr(M2-mm12-mm22)-4.*mm12*mm22)/(2.*sqrt(M2));
    p1   = Vec4D((M2-mm22+mm12)/(2.*sqrt(M2)),pabs*dir);
    p2   = Vec4D((M2+mm22-mm12)/(2.*sqrt(M2)),-pabs*dir);
    cluster->p_left->RescaleMomentum(p1);
    cluster->BoostBack(p2); cluster->BoostBack();
    Particle * h2 = new Particle(-1,cluster->p_right->FlavPair()->first,p2,'P');
    h2->SetNumber(0);
    p_hadrons->push_back(h2);
    testmom = cluster->Momentum()-cluster->p_left->Momentum()-p2;
  }
  else {
    msg.Error()<<"WARNING in Simple_Cluster_Fission::RearrangeMomenta: "<<std::endl
	       <<"   Seemingly only hadrons involved : "
	       <<cluster->p_left->FlavPair()->first<<" "<<cluster->p_right->FlavPair()->first
	       <<", continue with isotropic decay."<<std::endl;
    IsotropicDecay(cluster);
  }
}
 
void Simple_Cluster_Fission::IsotropicDecay(Cluster * cluster) 
{
  cluster->BoostInCMS();
  SelectHadronPair(cluster);
  Flavour flav1 = m_optiter->first.first, flav2 = m_optiter->first.second;
  double phi, cost, sint;
  Vec3D  dir;
  Particle * hadron, * part=NULL;
  if (cluster->GetParticle(0)->Info()=='L' || cluster->GetParticle(0)->Info()=='l') part = cluster->GetParticle(0);
  if (cluster->GetParticle(1)->Info()=='L' || cluster->GetParticle(1)->Info()=='l') part = cluster->GetParticle(1);
  if (part) {
    dir          = Vec3D(cluster->GetParticle(0)->Momentum()); dir = dir/dir.Abs();
    Poincare rot = Poincare(Vec4D(1.,dir),Vec4D(1.,0.,0.,1.));
    LeadingParticle(part,cost,phi);
    sint = sqrt(1-sqr(cost));
    dir          = Vec3D(sint*cos(phi),sint*sin(phi),cost);
    Vec4D help   = Vec4D(1.,dir);
    rot.RotateBack(help);
    dir          = Vec3D(help);
  }
  else {
    phi  = ran.Get()*2.*M_PI;
    cost = ran.Get()*2.-1.; sint = sqrt(1-sqr(cost));
    dir  = Vec3D(sint*cos(phi),sint*sin(phi),cost);
  }
  double m12  = sqr(flav1.Mass()), m22 = sqr(flav2.Mass());
  double pabs = sqrt(sqr(cluster->Mass2()-m12-m22)-4.*m12*m22)/(2.*cluster->Mass());
  double E1   = sqrt(sqr(pabs)+m12), E2 = sqrt(sqr(pabs)+m22);
  Vec4D  p1   = Vec4D(E1,pabs*dir), p2 = Vec4D(E2,-pabs*dir);
  cluster->BoostBack(p1); cluster->BoostBack(p2); cluster->BoostBack();
  hadron      = new Particle(-1,m_optiter->first.first,p1,'P');
  hadron->SetNumber(0);
  p_hadrons->push_back(hadron);
  hadron = new Particle(-1,m_optiter->first.second,p2,'P');
  hadron->SetNumber(0);
  p_hadrons->push_back(hadron);
}

void Simple_Cluster_Fission::LeadingParticle(ATOOLS::Particle * part,double & cost,double & phi)
{
  Flavour flav = part->Flav();
  double pref  = hadpars.GetConstituents()->Smearing(flav)*
    sqr(hadpars.GetConstituents()->Mass(flav)/part->Momentum()[0]);
  if (part->Info()=='l') pref *= 5.;
  cost  = 1.+pref*log(ran.Get());
  phi   = 2.*M_PI*ran.Get();
}

void Simple_Cluster_Fission::SelectHadronPair(Cluster * cluster) 
{
  m_options.clear();
  double      tot = 0, mass = cluster->Mass(), wt, m1, m2;
  FlavourPair pair;
  Allchanneliter clu = hadpars.GetAllChannels()->find(*(cluster->FlavPair()));
  if (clu==hadpars.GetAllChannels()->end()) {
    abort();
  }
  for (Channeliter ch=clu->second->begin();ch!=clu->second->end();ch++) {
    if (mass>((*ch)->m_mass)) {
      pair.first  = (*ch)->m_hadrons[0];
      pair.second = (*ch)->m_hadrons[1];
      
      m1          = (*ch)->m_hadrons[0].Mass();
      m2          = (*ch)->m_hadrons[1].Mass();
      wt          = sqr((*ch)->m_weight);
      wt         *= sqrt((sqr(mass)-sqr(m1+m2))*(sqr(mass)-sqr(m1-m2)));
      tot        += wt;
      m_options.insert(std::make_pair(pair,wt));
    }
  }
  if (m_options.empty()) {
    msg.Error()<<"Error in Simple_Cluster_Fission::SelectHadronPair."<<std::endl
	       <<"   Found no hadron decay option for cluster, will abort. Cluster : "<<std::endl
	       <<(*cluster)<<std::endl;
    abort();
  }
  tot *= ran.Get();
  for (m_optiter=m_options.begin();m_optiter!=m_options.end();m_optiter++) {
    tot -= m_optiter->second;
    if (tot<1.e-12) return;
  }
}

