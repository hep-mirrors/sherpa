#include "Interface_Tools.H"
#include "Random.H"

using namespace SHERPA;
using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


Interface_Tools::Interface_Tools(Tree ** _ini_trees,Tree * _fin_tree) :
  p_initrees(_ini_trees), p_fintree(_fin_tree) 
{ }

Interface_Tools::~Interface_Tools() { }

void Interface_Tools::InitializeIncoming(Blob * blob,double scale,
					 double th1,double th2)
{
  msg.Debugging()<<"In Interface_Tools::InitializeIncoming."<<std::endl;
  double E      = rpa.gen.Ecms();
  double x1     = (blob->InParton(0)->Momentum())[0]/E;
  double x2     = (blob->InParton(1)->Momentum())[0]/E;
  msg.Debugging()<<"Blob at E = "<<E<<" with x_{1,2} = "<<x1<<" ,  "<<x2<<std::endl
		 <<"Mom1 : "<<blob->InParton(0)->Momentum()<<std::endl
		 <<"Mom2 : "<<blob->InParton(1)->Momentum()<<std::endl;

  Knot * m1      = p_initrees[0]->NewKnot();
  *(m1->part)    = blob->InParton(0);
  m1->part->SetInfo('G');
  m1->part->SetStatus(1);
  m1->t          = -scale;
  m1->maxpt2     = scale;
  m1->costh      = -1.;
  m1->thcrit     = th1;
  m1->tout       = rpa.pshower.InitialQ02();
  m1->x          = x1;
  m1->E2         = sqr(x1*E);
  m1->stat       = 1;
  msg.Debugging()<<"First incoming parton : "<<m1<<std::endl;

  Knot * m2      = p_initrees[1]->NewKnot();
  *(m2->part)    = blob->InParton(1);
  m2->part->SetInfo('G');
  m2->part->SetStatus(1);
  m2->t          = -scale;
  m2->maxpt2     = scale;
  m2->costh      = -1.; 
  m2->thcrit     = th2;
  m2->tout       = rpa.pshower.InitialQ02();
  m2->x          = x2;
  m2->E2         = sqr(x2*E);
  m2->stat       = 1;
  msg.Debugging()<<"Second incoming parton : "<<m2<<std::endl;


  msg.Tracking()<<"Interface_Tools::InitializeIncoming :"<<std::endl
		<<"    "<<m1->kn_no<<" "<<m1->part->Flav()<<" "<<m1->E2<<" "<<m1->t<<" / "
		<<m2->kn_no<<" "<<m2->part->Flav()<<" "<<m2->E2<<" "<<m2->t<<std::endl;
}

void Interface_Tools::InitializeOutGoing(Blob * blob,double scale,
					 double th1,double th2)
{
  msg.Debugging()<<"In Interface_Tools::InitializeOutGoing :"<<scale<<std::endl;
  
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

  msg.Tracking()<<"Interface_Tools::InitializeOutGoing :"<<std::endl
		<<"    "<<d1->kn_no<<" "<<d1->part->Flav()<<" "<<d1->E2<<" "<<d1->t<<" / "
		<<d2->kn_no<<" "<<d2->part->Flav()<<" "<<d2->E2<<" "<<d2->t<<std::endl;
}

bool Interface_Tools::IsColourConnected(Parton * a, Parton * b) {
  msg.Debugging()<<"Check if "<<a->Flav()<<" and "<<b->Flav()<<" are connected."<<std::endl;
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
	msg.Debugging()<<"   partons "
		       <<a->Flav()<<"("<<a->Number()<<")  "
		       <<b->Flav()<<"("<<b->Number()<<")"<<std::endl
		       <<"   vecs= "<<avec<<" ,  "<<bvec<<std::endl;
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
	msg.Debugging()<<"   partons "
		       <<a->Flav()<<"("<<a->Number()<<")  "
		       <<b->Flav()<<"("<<b->Number()<<")"<<std::endl
		       <<"   vecs= "<<avec<<" ,  "<<bvec<<std::endl;
	angle = Max(angle,acos(avec*bvec/(avec.Abs()*bvec.Abs())));
      }
    }
  }
  msg.Debugging()<<"Angle = "<<angle<<std::endl;
  return angle;
}

