#include "Clusters_2_Hadrons.H"
#include "Message.H"
#include <iomanip>


using namespace AHADIC;
using namespace ATOOLS;

Clusters_2_Hadrons::Clusters_2_Hadrons() :
  m_Q(1.), p_newmasses(NULL), p_oldens2(NULL), p_ens(NULL),
  p_newflavs(NULL), p_newmoms(NULL)
{}

Clusters_2_Hadrons::~Clusters_2_Hadrons()
{
  if (p_newmasses) { delete [] p_newmasses; p_newmasses = NULL; }
  if (p_oldens2)   { delete [] p_oldens2;   p_oldens2 = NULL; }
  if (p_ens)       { delete [] p_ens;       p_ens = NULL; }
  if (p_newflavs)  { delete [] p_newflavs;  p_newflavs = NULL; }
  if (p_newmoms)   { delete [] p_newmoms;   p_newmoms   = NULL; }
}

int Clusters_2_Hadrons::Transition(Cluster * in)
{
  Allweightiter   all;
  Cluster * cluster=NULL;
  int             testit=0;
  double          mass,mass2,weight,maxweight,min;
  Flavour         hadron; 

  for (int i=0;i<2;i++) {
    if (i==0) cluster = in->p_left;
         else cluster = in->p_right;
    all     = hadpars.GetAllWeights()->find((*(cluster->FlavPair())));
    if (all==hadpars.GetAllWeights()->end()) {
      msg.Error()<<"Error in Clusters_2_Hadrons::Transition :"<<std::endl
		 <<"   Flavour combination {"<<cluster->FlavPair()->first<<"/"
		 <<cluster->FlavPair()->second<<"} "
		 <<" not found in list."<<std::endl;
      abort();
    }
    mass2   = cluster->Mass2();
    mass   = cluster->Mass();
    hadron = Flavour(kf::cluster);
    if (mass2 < 0. ||
	(mass < 1.01*all->second->begin()->first.Mass()) ||
	(mass < hadpars.MinimalMass(cluster->FlavPair())+m_Q) ||
	(mass < cluster->Mass(0)+cluster->Mass(1))) {
      maxweight = 0.;
      for (Weightiter test=all->second->begin();test!=all->second->end();test++) {
	weight  = 1./(sqr(sqr(mass)-sqr(test->first.Mass()))+0.001);
	weight *= sqr(test->second);
	if (weight>maxweight) { maxweight = weight; hadron = test->first; }
      }
      cluster->SetFlavour(hadron);
      mass    = hadron.Mass();
      testit += (i+1);
    }    
  }
  return testit;
}

void Clusters_2_Hadrons::Transition(Cluster_List * cl,Blob * blob)
{  
  if (cl->size()<2) return;
  Allweightiter   all;
  int             i=0;
  double          mass,weight,maxweight,min;
  Flavour         hadron; 

  Init(cl->size());
  for (Cluster_Iterator cit=cl->begin();cit!=cl->end();cit++,i++) {
    m_xmt        += p_newmasses[i] = sqrt((*cit)->Mass2());
    p_newmoms[i]  = (*cit)->Momentum();
    m_cms        += p_newmoms[i];
    p_oldens2[i]  = sqr(p_newmoms[i][0]);
  }
  ShiftBack(cl->size());

  i = 0;
  for (Cluster_Iterator cit=cl->begin();cit!=cl->end();cit++,i++) {
    all     = hadpars.GetAllWeights()->find((*((*cit)->FlavPair())));
    if (all==hadpars.GetAllWeights()->end()) {
      msg.Error()<<"Error in Clusters_2_Hadrons::Transition :"<<std::endl
		 <<"   Flavour combination {"<<(*cit)->FlavPair()->first<<"/"
		 <<(*cit)->FlavPair()->second<<"} "
		 <<" not found in list."<<std::endl;
      abort();
    }
    mass   = sqrt((*cit)->Mass2());
    hadron = Flavour(kf::cluster);
    if ((mass < 1.01*all->second->begin()->first.Mass()) ||
	(mass < hadpars.MinimalMass((*cit)->FlavPair())+m_Q) ||
	(mass < (*cit)->Mass(0)+(*cit)->Mass(1))) {
      maxweight = 0.;
      for (Weightiter test=all->second->begin();test!=all->second->end();test++) {
	weight  = 1./(sqr(sqr(mass)-sqr(test->first.Mass()))+0.001);
	weight *= sqr(test->second);
	if (weight>maxweight) { maxweight = weight; hadron = test->first; }
      }
      mass = hadron.Mass();
    }
    p_newflavs[i] = hadron;
    m_xmt        += p_newmasses[i] = mass;
    p_oldens2[i]  = sqr(p_newmoms[i][0]);
  }
  Shift(cl->size());

  i = 0;
  Cluster_Iterator cit = cl->begin();
  do {
    if (p_newflavs[i].IsHadron()) {
      AddToBlob(blob,i);
      if (*cit) delete (*cit);
      cit = cl->erase(cit);
    } 
    else {
      (*cit)->RescaleMomentum(p_newmoms[i]);
      cit++;
    }
    i++;
  } while (cit!=cl->end());
}


void Clusters_2_Hadrons::AddToBlob(ATOOLS::Blob * blob,const int pos)
{
  Particle * part = new Particle(-1,p_newflavs[pos],p_newmoms[pos],'P');
  part->SetNumber(0);
  blob->AddToOutParticles(part);
}


void Clusters_2_Hadrons::ShiftBack(size_t number)
{
  double ET  = sqrt(m_cms.Abs2()); 
  double x   = sqrt(1.-sqr(m_xmt/ET));
  double acc = ET*1.e-14;

  double f0,g0,x2;
  for (int j=0;j<10;j++) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=0;i<number;i++) {
      p_ens[i] = sqrt((p_oldens2[i]-sqr(p_newmasses[i]))*x2); //);
      f0      += p_ens[i];
      g0      += (p_oldens2[i]-sqr(p_newmasses[i]))/p_ens[i];
    }
    if (dabs(f0)<acc) break; 
    x -= f0/(g0*x);  
  }
  
  for (int j=0;j<number;j++) p_newmoms[j] = Vec4D(p_ens[j],x*Vec3D(p_newmoms[j]));
  m_xmt = 0.;
}

void Clusters_2_Hadrons::Shift(size_t number)
{
  double ET  = sqrt(m_cms.Abs2()); 
  double x   = sqrt(1.-sqr(m_xmt/ET));
  double acc = ET*1.e-14;

  double f0,g0,x2;
  for (int j=0;j<10;j++) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=0;i<number;i++) {
      p_ens[i] = sqrt(sqr(p_newmasses[i])+x2*p_oldens2[i]);
      f0      += p_ens[i];
      g0      += p_oldens2[i]/p_ens[i];
    }
    if (dabs(f0)<acc) break; 
    x -= f0/(x*g0);  
  }
  // Construct Momenta
  
  for (int j=0;j<number;j++) p_newmoms[j] = Vec4D(p_ens[j],x*Vec3D(p_newmoms[j]));
}

void Clusters_2_Hadrons::Init(int number)
{
  if (p_newmasses) delete [] p_newmasses; 
  if (p_oldens2)   delete [] p_oldens2; 
  if (p_ens)       delete [] p_ens;     
  if (p_newflavs)  delete [] p_newflavs;
  if (p_newmoms)   delete [] p_newmoms; 

  p_newmasses = new double[number],
  p_oldens2   = new double[number],
  p_ens       = new double[number];
  p_newflavs  = new Flavour[number];
  p_newmoms   = new Vec4D[number];

  m_xmt       = 0.; 
  m_cms       = Vec4D(0.,0.,0.,0.);
}



