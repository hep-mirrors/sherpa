#include "Interface_Tools.H"

#include "Random.H"

using namespace SHERPA;
using namespace APACIC;
using namespace ATOOLS;

Interface_Tools::Interface_Tools(Tree **ini_trees,Tree *fin_tree):
  p_initrees(ini_trees), 
  p_fintree(fin_tree), 
  p_dummy(new Particle(0,Flavour(kf::none),ATOOLS::Vec4D())) {}
	  
Interface_Tools::~Interface_Tools() 
{ 
  delete p_dummy;
}

void Interface_Tools::InitializeIncoming(Blob *blob,double scale,double E,
					 double th1,double th2,double x1,double x2)
{
  p_initrees[0]->Reset();  
  Knot * m1      = p_initrees[0]->NewKnot();
  *(m1->part)    = blob->InParticle(0);
  m1->part->SetInfo('G');
  m1->part->SetStatus(1);
  m1->t          = -scale;
  m1->maxpt2     = scale;
  m1->costh      = -1.;
  m1->thcrit     = th1;
  m1->tout       = sqr((m1->part->Flav()).PSMass());
  m1->x          = x1;
  m1->E2         = sqr(x1*E);
  m1->stat       = 1;
  m1->part->SetDecayBlob(blob);
  p_initrees[1]->Reset();  
  Knot * m2      = p_initrees[1]->NewKnot();
  *(m2->part)    = blob->InParticle(1);
  m2->part->SetInfo('G');
  m2->part->SetStatus(1);
  m2->t          = -scale;
  m2->maxpt2     = scale;
  m2->costh      = -1.; 
  m2->thcrit     = th2;
  m2->tout       = sqr((m2->part->Flav()).PSMass());
  m2->x          = x2;
  m2->E2         = sqr(x2*E);
  m2->stat       = 1;
  m2->part->SetDecayBlob(blob);
}

void Interface_Tools::InitializeOutGoing(Blob *blob,double scale,double E,
					 double th1,double th2)
{
  blob->SetCMS();
  blob->BoostInCMS();
  ATOOLS::Particle *part1=blob->OutParticle(0);
  ATOOLS::Particle *part2=blob->OutParticle(1);
  p_fintree->Reset();
  p_dummy->SetMomentum(part1->Momentum()+part2->Momentum());
  Knot * dummy   = p_fintree->NewKnot(p_dummy);
  dummy->part->SetInfo('M');
  dummy->part->SetStatus(2);
  dummy->t       = scale;
  dummy->maxpt2  = scale;
  dummy->costh   = -1;
  dummy->thcrit  = M_PI;
  dummy->stat    = 0;
  dummy->E2      = sqr(dummy->part->Momentum()[0]);
  Knot * d1      = p_fintree->NewKnot(part1);
  d1->part->SetInfo('H');
  d1->part->SetStatus(1);
  d1->t          = scale;
  d1->maxpt2     = scale;
  d1->costh      = -1.; 
  d1->thcrit     = 1.e50;
  d1->tout       = sqr((d1->part->Flav()).PSMass());
  d1->E2         = sqr(d1->part->Momentum()[0]);
  d1->stat       = 3;
  d1->part->SetProductionBlob(blob);
  Knot * d2      = p_fintree->NewKnot(part2);
  d2->part->SetInfo('H');
  d2->part->SetStatus(1);
  d2->t          = scale;
  d2->maxpt2     = scale;
  d2->costh      = -1.; 
  d2->thcrit     = 1.e50;
  d2->tout       = sqr((d2->part->Flav()).PSMass());
  d2->E2         = sqr(d2->part->Momentum()[0]);
  d2->stat       = 3;
  d2->part->SetProductionBlob(blob);
  dummy->E2      = sqr(sqrt(d1->E2)+sqrt(d2->E2));
  dummy->z       = sqrt(d1->E2/dummy->E2);
  d1->prev       = dummy;
  dummy->left    = d1;
  d2->prev       = dummy;
  dummy->right   = d2;
  blob->BoostInLab();
}

bool Interface_Tools::IsColourConnected(Particle * a, Particle * b) {
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				 (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				 (a->GetFlow(2)==b->GetFlow(1)))  )    );
}

double Interface_Tools::ColourAngle(Particle * a,Blob * blob) {
  double angle = M_PI;
  if (!((a->Flav()).Strong())) return angle;

  Particle * b;
  angle = 0.;
  Vec3D avec,bvec;
  for (int i=0;i<blob->NInP();i++) {
    b = blob->InParticle(i);
    if (b != a) {
      if (IsColourConnected(a,b)) {
	avec = Vec3D(a->Momentum());
	bvec = Vec3D(b->Momentum());
	angle = Max(angle,acos(avec*bvec/(avec.Abs()*bvec.Abs())));
      }
    }
  }
  for (int i=0;i<blob->NOutP();i++) {
    b = blob->OutParticle(i);
    if (b != a) {
      if (IsColourConnected(a,b)) {
	avec = Vec3D(a->Momentum());
	bvec = Vec3D(b->Momentum());
	angle = Max(angle,acos(avec*bvec/(avec.Abs()*bvec.Abs())));
      }
    }
  }
  return angle;
}

