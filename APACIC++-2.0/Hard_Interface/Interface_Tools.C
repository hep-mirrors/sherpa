#include "Interface_Tools.H"
#include "Random.H"

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMEGIC;


Interface_Tools::Interface_Tools(Tree ** _ini_trees,Tree * _fin_tree) :
  ini_trees(_ini_trees), fin_tree(_fin_tree)
{ }

Interface_Tools::~Interface_Tools() 
{ }

void Interface_Tools::InitializeIncoming(Blob * blob,double scale,
					 double th1,double th2)
{
  msg.Debugging()<<"In Interface_Tools::InitializeIncoming."<<std::endl;
  double E      = rpa.gen.Ecms();
  double x1     = (blob->InParton(0)->momentum())[0]/E;
  double x2     = (blob->InParton(1)->momentum())[0]/E;
  msg.Debugging()<<"Blob at E = "<<E<<" with x_{1,2} = "<<x1<<" ,  "<<x2<<std::endl
		 <<"Mom1 : "<<blob->InParton(0)->momentum()<<std::endl
		 <<"Mom2 : "<<blob->InParton(1)->momentum()<<std::endl;

  Knot * m1      = ini_trees[0]->NewKnot();
  *(m1->part)    = blob->InParton(0);
  m1->part->set_info('G');
  m1->part->set_status(1);
  m1->t          = -scale;
  m1->maxpt2     = scale;
  m1->costh      = -1.;
  m1->thcrit     = th1;
  m1->tout       = rpa.pshower.InitialQ02();
  m1->x          = x1;
  m1->E2         = sqr(x1*E);
  m1->stat       = 1;
  msg.Debugging()<<"First incoming parton : "<<m1<<std::endl;

  Knot * m2      = ini_trees[1]->NewKnot();
  *(m2->part)    = blob->InParton(1);
  m2->part->set_info('G');
  m2->part->set_status(1);
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
		<<"    "<<m1->kn_no<<" "<<m1->part->flav()<<" "<<m1->E2<<" "<<m1->t<<" / "
		<<m2->kn_no<<" "<<m2->part->flav()<<" "<<m2->E2<<" "<<m2->t<<std::endl;
}

void Interface_Tools::InitializeOutGoing(Blob * blob,double scale,
					 double th1,double th2)
{
  msg.Debugging()<<"In Interface_Tools::InitializeOutGoing :"<<scale<<std::endl;
  
  Knot * dummy   = fin_tree->NewKnot();
  *(dummy->part) = Parton(0,Flavour(kf::none),blob->CMS());
  dummy->part->set_info('M');
  dummy->part->set_status(2);
  dummy->t       = scale;
  dummy->maxpt2  = scale;
  dummy->costh   = -1;
  dummy->thcrit  = M_PI;
  dummy->stat    = 1;
    
  Knot * d1      = fin_tree->NewKnot();
  *(d1->part)    = blob->OutParton(0);
  d1->part->set_info('H');
  d1->part->set_status(1);
  d1->t          = scale;
  d1->maxpt2     = scale;
  d1->costh      = -1.; 
  d1->thcrit     = th1;
  d1->tout       = sqr((d1->part->flav()).PSmass());
  d1->E2         = sqr(d1->part->momentum()[0]);
  d1->stat       = 1;

  Knot * d2      = fin_tree->NewKnot();
  *(d2->part)    = blob->OutParton(1);
  d2->part->set_info('H');
  d2->part->set_status(1);
  d2->t          = scale;
  d2->maxpt2     = scale;
  d2->costh      = -1.; 
  d2->thcrit     = th2;
  d2->tout       = sqr((d2->part->flav()).PSmass());
  d2->E2         = sqr(d2->part->momentum()[0]);
  d2->stat       = 1;
  
  dummy->E2      = sqr(sqrt(d1->E2)+sqrt(d2->E2));
  dummy->z       = sqrt(d1->E2/dummy->E2);
  
  d1->prev       = dummy;
  dummy->left    = d1;
  d2->prev       = dummy;
  dummy->right   = d2;

  msg.Tracking()<<"Interface_Tools::InitializeOutGoing :"<<std::endl
		<<"    "<<d1->kn_no<<" "<<d1->part->flav()<<" "<<d1->E2<<" "<<d1->t<<" / "
		<<d2->kn_no<<" "<<d2->part->flav()<<" "<<d2->E2<<" "<<d2->t<<std::endl;
}

bool Interface_Tools::IsColourConnected(Parton * a, Parton * b) {
  msg.Debugging()<<"Check if "<<a->flav()<<" and "<<b->flav()<<" are connected."<<std::endl;
  return (( (a->flow(1)!=0) && ( (a->flow(1)==b->flow(1)) || 
				 (a->flow(1)==b->flow(2)))  ) ||
	  ( (a->flow(2)!=0) && ( (a->flow(2)==b->flow(2)) ||
				 (a->flow(2)==b->flow(1)))  )    );
}

double Interface_Tools::ColourAngle(Parton * a,Blob * blob) {
  double angle = M_PI;
  if (!((a->flav()).strong())) return angle;

  Parton * b;
  angle = 0.;
  vec3d avec,bvec;
  for (int i=0;i<blob->NInP();i++) {
    b = blob->InParton(i);
    if (b != a) {
      if (IsColourConnected(a,b)) {
	avec = vec3d(a->momentum());
	bvec = vec3d(b->momentum());
	msg.Debugging()<<"   partons "
		       <<a->flav()<<"("<<a->Get_Numb()<<")  "
		       <<b->flav()<<"("<<b->Get_Numb()<<")"<<std::endl
		       <<"   vecs= "<<avec<<" ,  "<<bvec<<std::endl;
	angle = Max(angle,acos(avec*bvec/(avec.abs()*bvec.abs())));
      }
    }
  }
  for (int i=0;i<blob->NOutP();i++) {
    b = blob->OutParton(i);
    if (b != a) {
      if (IsColourConnected(a,b)) {
	avec = vec3d(a->momentum());
	bvec = vec3d(b->momentum());
	msg.Debugging()<<"   partons "
		       <<a->flav()<<"("<<a->Get_Numb()<<")  "
		       <<b->flav()<<"("<<b->Get_Numb()<<")"<<std::endl
		       <<"   vecs= "<<avec<<" ,  "<<bvec<<std::endl;
	angle = Max(angle,acos(avec*bvec/(avec.abs()*bvec.abs())));
      }
    }
  }
  msg.Debugging()<<"Angle = "<<angle<<std::endl;
  return angle;
}

