#include "Interface_Tools.H"
#include "Random.H"

using namespace SHERPA;
using namespace APACIC;
using namespace ATOOLS;


Interface_Tools::Interface_Tools(Tree ** _ini_trees,Tree * _fin_tree) :
  p_initrees(_ini_trees), p_fintree(_fin_tree) 
{ }

Interface_Tools::~Interface_Tools() { }

void Interface_Tools::InitializeIncoming(Blob * blob,double scale,double E,
					 double th1,double th2,double x1,double x2)
{
  Knot * m1      = p_initrees[0]->NewKnot();
  *(m1->part)    = blob->InParton(0);
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

  Knot * m2      = p_initrees[1]->NewKnot();
  *(m2->part)    = blob->InParton(1);
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
}

void Interface_Tools::InitializeOutGoing(Blob * blob,double scale,double E,
					 double th1,double th2)
{
  Knot * dummy   = p_fintree->NewKnot(new Parton(0,Flavour(kf::none),blob->CMS()));
  dummy->part->SetInfo('M');
  dummy->part->SetStatus(2);
  dummy->t       = scale;
  dummy->maxpt2  = scale;
  dummy->costh   = -1;
  dummy->thcrit  = M_PI;
  dummy->stat    = 1;
    
  Knot * d1      = p_fintree->NewKnot(blob->OutParton(0));
  d1->part->SetInfo('H');
  d1->part->SetStatus(1);
  d1->t          = scale;
  d1->maxpt2     = scale;
  d1->costh      = -1.; 
  d1->thcrit     = 1.e50;
  d1->tout       = sqr((d1->part->Flav()).PSMass());
  d1->E2         = sqr(d1->part->Momentum()[0]);
  d1->stat       = 1;

  Knot * d2      = p_fintree->NewKnot(blob->OutParton(1));
  d2->part->SetInfo('H');
  d2->part->SetStatus(1);
  d2->t          = scale;
  d2->maxpt2     = scale;
  d2->costh      = -1.; 
  d2->thcrit     = 1.e50;
  d2->tout       = sqr((d2->part->Flav()).PSMass());
  d2->E2         = sqr(d2->part->Momentum()[0]);
  d2->stat       = 1;
  
  dummy->E2      = sqr(sqrt(d1->E2)+sqrt(d2->E2));
  dummy->z       = sqrt(d1->E2/dummy->E2);
  
  d1->prev       = dummy;
  dummy->left    = d1;
  d2->prev       = dummy;
  dummy->right   = d2;
}

bool Interface_Tools::IsColourConnected(Parton * a, Parton * b) {
  return (( (a->GetFlow(1)!=0) && ( (a->GetFlow(1)==b->GetFlow(1)) || 
				 (a->GetFlow(1)==b->GetFlow(2)))  ) ||
	  ( (a->GetFlow(2)!=0) && ( (a->GetFlow(2)==b->GetFlow(2)) ||
				 (a->GetFlow(2)==b->GetFlow(1)))  )    );
}

double Interface_Tools::ColourAngle(Parton * a,Blob * blob) {
  double angle = M_PI;
  if (!((a->Flav()).Strong())) return angle;

  Parton * b;
  angle = 0.;
  Vec3D avec,bvec;
  for (int i=0;i<blob->NInP();i++) {
    b = blob->InParton(i);
    if (b != a) {
      if (IsColourConnected(a,b)) {
	avec = Vec3D(a->Momentum());
	bvec = Vec3D(b->Momentum());
	angle = Max(angle,acos(avec*bvec/(avec.Abs()*bvec.Abs())));
      }
    }
  }
  for (int i=0;i<blob->NOutP();i++) {
    b = blob->OutParton(i);
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

