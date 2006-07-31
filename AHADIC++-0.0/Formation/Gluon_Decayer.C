#include "Gluon_Decayer.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;


Gluon_Decayer::Gluon_Decayer()
{ 
  Flavour flav = Flavour(kf::gluon);
  double gluonmass = hadpars.GetConstituents()->Mass(flav);
  for (FlavCCMap_Iterator fdit=hadpars.GetConstituents()->CCMap.begin();
       fdit!=hadpars.GetConstituents()->CCMap.end();fdit++) {
    if (2.*fdit->second->Mass()<gluonmass && !fdit->first.IsAnti()) {
      m_options.insert(std::make_pair(fdit->first,new DecaySpecs()));
      msg.Tracking()<<"Insert option : g->"<<fdit->first<<" "<<fdit->first.Bar()<<std::endl;
    }
  }
  msg.Tracking()<<"------------- END OF GLUON_DECAYER --------------"<<std::endl;
}

Gluon_Decayer::~Gluon_Decayer() {
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    if (fdit->second!=NULL) { delete fdit->second; fdit->second=NULL; }
  }
  m_options.clear();
}

bool Gluon_Decayer::DecayList(Part_List * plin)
{
  if (plin==NULL || plin->empty()) return true;

  //Vec4D mom=Vec4D(0.,0.,0.,0.);
  //for (Part_Iterator pit=plin->begin();pit!=plin->end();pit++) mom+=(*pit)->Momentum();
  //std::cout<<"Before : "<<plin->size()<<" with "<<mom<<"."<<std::endl;
  //if (plin->size()==2) {
  //  for (Part_Iterator pit=plin->begin();pit!=plin->end();pit++) 
  //    std::cout<<"   "<<(*pit)->Flav()<<" : "<<(*pit)->Momentum()<<std::endl;
  //}

  if (!Shift(plin)) return false;

  //mom=Vec4D(0.,0.,0.,0.);
  //for (Part_Iterator pit=plin->begin();pit!=plin->end();pit++) mom+=(*pit)->Momentum();
  //std::cout<<"Shift  : "<<plin->size()<<" with "<<mom<<"."<<std::endl;

  bool success(true);
  Particle * part, * part1, * part2;
  Part_Iterator plit=plin->begin();
  do {
    part = (*plit);
    if (part && part->Info()=='L') {
      part1 = part2 = NULL;
      if (!DecayIt(part,part1,part2)) success=false;
      else {
	if (part1) {
	  delete part;
	  plit = plin->erase(plit);
	  plit = plin->insert(plit,part1);
	  plit = plin->insert(plit,part2);
	  ++plit;++plit;
	}
	else ++plit;
      }
    } 
    else ++plit;
  } while (plit!=plin->end());

  //mom=Vec4D(0.,0.,0.,0.);
  //for (Part_Iterator pit=plin->begin();pit!=plin->end();pit++) mom+=(*pit)->Momentum();
  //std::cout<<"After  : "<<plin->size()<<" with "<<mom<<"."<<std::endl;

  return success;
}

bool Gluon_Decayer::Shift(Part_List * pl)
{
  size_t number(pl->size());
  bool val(true);
  if (number<2) return val; 
  Vec4D  * momenta = new Vec4D[number];
  double * masses  = new double[number];
  int k(0);
  Flavour flav;
  for (Part_Iterator pit=pl->begin();pit!=pl->end();++pit,++k) {
    flav       = (*pit)->Flav();
    /*
      if (flav==Flavour(kf::c)  || flav==Flavour(kf::c).Bar() ||
      flav==Flavour(kf::b)  || flav==Flavour(kf::b).Bar()) {
      bool bar = flav.IsAnti();
      std::cout<<"Change heavy flavour : "<<flav<<" -> ";
      flav     = Flavour(int(1+3*ran.Get()));
      if (bar) flav = flav.Bar();
      (*pit)->SetFlav(flav);
      std::cout<<flav<<std::endl;
      }
    */
    momenta[k] = (*pit)->Momentum();
    masses[k]  = hadpars.GetConstituents()->Mass(flav);
    //std::cout<<"Gluon_Decayer : "<<k<<" : "<<flav<<" "<<masses[k]<<std::endl;
  }
  if (!hadpars.AdjustMomenta(number,momenta,masses)) val=false;

  k = 0;
  for (Part_Iterator pit=pl->begin();pit!=pl->end();++pit,++k) 
    (*pit)->SetMomentum(momenta[k]);

  delete momenta;
  delete masses;
  return val;
}

bool Gluon_Decayer::DecayIt(ATOOLS::Particle * part,
			    ATOOLS::Particle *& part1,ATOOLS::Particle *& part2)
{
  if (!part->Flav().IsGluon()) return true;
  SelectDecay(part->Momentum());
  BuildKinematics(part->Momentum());
  part1 = new Particle(0,m_flav,m_p1vec,'l');
  part1->SetNumber(0);
  part1->SetFlow(1,part->GetFlow(1));
  part1->SetStatus(part_status::active);

  part2 = new Particle(0,m_flav.Bar(),m_p2vec,'l');
  part2->SetNumber(0);
  part2->SetFlow(2,part->GetFlow(2));
  part2->SetStatus(part_status::active);
  return true;
}

void Gluon_Decayer::SelectDecay(const Vec4D & mom) {
  std::vector<double> weights;
  double  disc,sum=0.,mass;
  Flavour flav;
  FDIter  fdit;
  for (fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    flav   = fdit->first;
    mass   = hadpars.GetConstituents()->Mass(flav);
    disc   = Vec3D(mom).Abs()/mom[0]*
      sqrt(sqr(mom.Abs2()-2.*sqr(mass))-4.*pow(mass,4))/mom.Abs2();
    fdit->second->zmin   = (1-disc)/2.;
    fdit->second->zmax   = (1+disc)/2.;
    sum   += fdit->second->weight = ZWeight(fdit->second->zmin,
					    fdit->second->zmax);
    weights.push_back(fdit->second->weight);
  }

  disc = sum*ran.Get();
  for (fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    disc -= fdit->second->weight;
    if (disc<0.) break;
  }
  m_flav = fdit->first;
  m_mass = hadpars.GetConstituents()->Mass(m_flav);
  SelectZ(fdit->second);
  SelectPhi();
}


void Gluon_Decayer::BuildKinematics(const Vec4D & mom) 
{
  double E  = mom[0],      p  = Vec3D(mom).Abs();
  double t1 = sqr(m_mass), E1 = m_z*E,      p1 = sqrt(E1*E1-t1);
  double t2 = t1,          E2 = (1.-m_z)*E, p2 = sqrt(E2*E2-t2);
    
  double cth1 = (p*p-p2*p2+p1*p1)/(2.*p*p1), sth1 = sqrt(1.-sqr(cth1));
  double cth2 = (p*p+p2*p2-p1*p1)/(2.*p*p2), sth2 = sqrt(1.-sqr(cth2));

  Vec3D nm = Vec3D(mom)/Vec3D(mom).Abs();
  Vec3D n1 = cross(Vec3D(0.,0.,1.),nm);
  if (n1.Abs()<1.e-5) {
    n1     = cross(Vec3D(0.,1.,0.),nm);
  }
  n1       = n1/n1.Abs();
  Vec3D es = cos(m_phi)*n1 + sin(m_phi)*cross(nm,n1);

  m_p1vec  = Vec4D(E1,p1*(cth1*nm - sth1*es));
  m_p2vec  = Vec4D(E2,p2*(cth2*nm + sth2*es));
}


double Gluon_Decayer::ZWeight(const double zmin,const double zmax) {
  return (zmax-zmin)-sqr(zmax-zmin)+2.*pow(zmax-zmin,3)/3.;
}

void Gluon_Decayer::SelectZ(const DecaySpecs * ds) 
{
  for (;;) {
    m_z = ds->zmin+ran.Get()*(ds->zmax-ds->zmin);
    if (sqr(m_z)+sqr(1-m_z)>ran.Get()) return;
  }
}

void Gluon_Decayer::SelectPhi()
{
  for (;;) {
    m_phi = ran.Get()*2.*M_PI;
    return;
  }
}
