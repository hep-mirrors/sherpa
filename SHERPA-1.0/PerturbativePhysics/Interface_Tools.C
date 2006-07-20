#include "Interface_Tools.H"

#include "Particle.H"
#include "Blob.H"
#include "Tree.H"
#include "Random.H"

using namespace SHERPA;
using namespace APACIC;
using namespace ATOOLS;

Interface_Tools::Interface_Tools(Tree **ini_trees,Tree *fin_tree):
  p_initrees(ini_trees), 
  p_fintree(fin_tree) {}
	  
Interface_Tools::~Interface_Tools() 
{ 
}

void Interface_Tools::InitializeIncoming(const Blob *blob,const double &E)
{
  const ATOOLS::Particle *part1=blob->ConstOutParticle(0);
  const ATOOLS::Particle *part2=blob->ConstOutParticle(1);
  double scale=(part1->Momentum()+part2->Momentum()).Abs2();
  Knot *m1=p_initrees[0]->NewKnot();
  *(m1->part)=*blob->ConstInParticle(0);
  m1->part->SetInfo('G');
  m1->part->SetStatus(1);
  m1->t=-scale;
  m1->costh=-1.;
  m1->thcrit=Angle(blob->ConstInParticle(0),blob);
  m1->tout=sqr(m1->part->Flav().PSMass());
  m1->x=blob->ConstInParticle(0)->Momentum()[0]/E;
  m1->E2=sqr(m1->x*E);
  m1->stat=1;
  m1->part->SetDecayBlob((ATOOLS::Blob*)blob);
  Knot *m2=p_initrees[1]->NewKnot();
  *(m2->part)=*blob->ConstInParticle(1);
  m2->part->SetInfo('G');
  m2->part->SetStatus(1);
  m2->t=-scale;
  m2->costh=-1.; 
  m2->thcrit=Angle(blob->ConstInParticle(0),blob);
  m2->tout=sqr(m2->part->Flav().PSMass());
  m2->x=blob->ConstInParticle(1)->Momentum()[0]/E;
  m2->E2=sqr(m2->x*E);
  m2->stat=1;
  m2->part->SetDecayBlob((ATOOLS::Blob*)blob);
  m_inipt2=part1->Momentum().PPerp2();
  m1->maxpt2=m1->pt2lcm=m_inipt2;
  m2->maxpt2=m2->pt2lcm=m_inipt2;
}

void Interface_Tools::InitializeOutGoing(Blob *blob,const double &E)
{
  ATOOLS::Particle *part1=blob->OutParticle(0);
  ATOOLS::Particle *part2=blob->OutParticle(1);
  Knot *dummy=p_fintree->NewKnot();
  dummy->part->SetMomentum(part1->Momentum()+part2->Momentum());
  Knot *d1=p_fintree->NewKnot(part1);
  Knot *d2=p_fintree->NewKnot(part2);
  blob->SetCMS();
  blob->BoostInCMS();
  dummy->part->SetInfo('f');
  dummy->part->SetStatus(2);
  dummy->t=dummy->maxpt2=dummy->tout=dummy->part->Momentum().Abs2();
  dummy->costh=-1;
  dummy->thcrit=M_PI;
  dummy->stat=0;
  dummy->E2=sqr(dummy->part->Momentum()[0]);
  dummy->zs=dummy->z=part1->Momentum()[0]/dummy->part->Momentum()[0];
  dummy->didkin=true;
  d1->part->SetInfo('H');
  d1->part->SetStatus(1);
  d1->t=dummy->t;
  d1->costh=-1.; 
  d1->thcrit=Angle(blob->InParticle(0),blob);
  d1->tout=sqr(part1->Flav().PSMass());
  d1->E2=sqr(part1->Momentum()[0]);
  d1->stat=3;
  d1->part->SetProductionBlob(blob);
  d1->didkin=true;
  d2->part->SetInfo('H');
  d2->part->SetStatus(1);
  d2->t=dummy->t;
  d2->costh=-1.; 
  d2->thcrit=Angle(blob->InParticle(1),blob);
  d2->tout=sqr(part2->Flav().PSMass());
  d2->E2=sqr(part2->Momentum()[0]);
  d2->stat=3;
  d2->part->SetProductionBlob(blob);
  d2->didkin=true;
  d1->prev=d2->prev=dummy;
  dummy->left=d1;
  dummy->right=d2;
  blob->BoostInLab();
  m_finpt2=part1->Momentum().PPerp2();
  d1->maxpt2=d1->pt2lcm=m_finpt2;
  d2->maxpt2=d2->pt2lcm=m_finpt2;
}

bool Interface_Tools::Connected(const Particle *a,const Particle *b) 
{
  return (a->GetFlow(1)!=0 && 
	  (a->GetFlow(1)==b->GetFlow(1) || a->GetFlow(1)==b->GetFlow(2))) ||
    (a->GetFlow(2)!=0 && 
     (a->GetFlow(2)==b->GetFlow(2) || a->GetFlow(2)==b->GetFlow(1)));
}

double Interface_Tools::Angle(const Particle *p1,const Blob *blob) 
{
  if (!p1->Flav().Strong()) return M_PI;
  double angle=0.0;
  for (int i=0;i<blob->NInP();++i) {
    const Particle *p2=blob->ConstInParticle(i);
    if (p2!=p1 && Connected(p2,p1)) 
      angle=ATOOLS::Max(angle,p2->Momentum().Theta(p1->Momentum()));
  }
  for (int i=0;i<blob->NOutP();i++) {
    const Particle *p2=blob->ConstOutParticle(i);
    if (p2!=p1 && Connected(p2,p1)) 
      angle=ATOOLS::Max(angle,p2->Momentum().Theta(p1->Momentum()));
  }
  return angle;
}

void Interface_Tools::JetVetoPt2(double &inipt2,double &finpt2)
{
  inipt2=m_inipt2;
  finpt2=m_finpt2;
}
