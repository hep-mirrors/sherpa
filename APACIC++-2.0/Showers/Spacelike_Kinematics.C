#include "Spacelike_Kinematics.H"
#include "Run_Parameter.H"
#include <iomanip.h>

using std::endl;

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;


void Spacelike_Kinematics::InitKinematics(Knot *& k1, Knot *& k2) {
  if ((!k1) || (!k2)) {
    msg.Error()<<"ERROR in Spacelike_Kinematics::InitKinematics : No knots."<<endl;
    return;
  }  
  rot           = Poincare(ZVEC,k1->part->Momentum());
  Vec4D cms     = k1->part->Momentum() + k2->part->Momentum();
  double sprime = cms.Abs2();
  double E1     = (sprime + k1->t - k2->t)/sqrt(4.*sprime);
  double E2     = (sprime - k1->t + k2->t)/sqrt(4.*sprime);
  double pz     = sqrt((sqr(sprime - k1->t - k2->t)-4.*k1->t*k2->t)/(4.*sprime));
  Vec4D v1(E1,0.,0.,pz);
  Vec4D v2(E2,0.,0.,-pz);
  k1->part->SetMomentum(v1);
  k2->part->SetMomentum(v2);  

  msg.Debugging()<<"Spacelike_Kinematics::InitKinematics :"<<endl
		 <<"   Vec1 : "<<v1<<" : "<<v1.Abs2()<<" / "<<k1->t<<endl
		 <<"   Vec2 : "<<v2<<" : "<<v2.Abs2()<<" / "<<k2->t<<endl
		 <<"   S    : "<<(v1+v2).Abs2()<<" / "<<sprime<<endl;
}

bool Spacelike_Kinematics::DoKinematics(Tree ** trees,Knot * active, Knot * partner,int leg) {
  if (!active->prev) {
    msg.Error()<<"Error Spacelike_Kinematics::DoKinematics : "
	       <<"     No mother for active knot, no kinematics to be constructed"<<endl;
    return 0;
  }

  BoostInCMS(trees,active, partner);

  Vec4D cms      = active->part->Momentum()+partner->part->Momentum();
  double sprime  = cms.Abs2();

  Knot * mother  = active->prev;
  Knot * sister  = mother->left;

  double s1      = sprime           - (partner->t) - (active->t);    
  double s3      = sprime/active->z - (partner->t) - (mother->t);
  double np1     = sqrt(s1*s1-4.*(partner->t)*(active->t));
  double np3     = sqrt(s3*s3-4.*(partner->t)*(mother->t));
  double maxt_d2 = CalculateMaxT(active,partner);

  if (maxt_d2 < sister->t) return 0;
  
  double E_mo      = 1./(2.*sqrt(sprime)) *
    (sprime/active->z - partner->t + active->t - sister->t);
  double pz_mo     = 1./(2.*active->part->Momentum()[3]) * 
    (s3 - 2. * partner->part->E() * E_mo);
  double pt_mo     = sqrt((maxt_d2-sister->t)/(np1*np1) *
    (0.5*(s1*s3 + np1*np3) + partner->t * (sister->t - active->t - mother->t)));
 
  msg.Debugging()<<"Spacelike_Kinematics::DoKinematics :"<<endl
		 <<"   CMS, s' = "<<cms<<", "<<sprime<<endl
		 <<"   timelike daughter's t : "<<maxt_d2<<"  >?>  "<<sister->t<<endl;

  double cph=cos(active->phi), sph=sin(active->phi);
  
  Vec4D v_mo(E_mo,sph*pt_mo,cph*pt_mo,pz_mo);
  Vec4D v_si = v_mo + (-1.)*active->part->Momentum();

  bool error = 0;
  // Checks:
  error = CheckVector(active->part->Momentum()) || CheckVector(partner->part->Momentum()) ||
    CheckVector(v_mo) || CheckVector(v_si);

  if (error) {
    msg.Error()<<"Error in Spacelike_Kinematics::DoKinematics : Bad vectors. "<<endl
	       <<"  Act : "<<active->part->Momentum()<<" : "<<active->part->Momentum().Abs2()
	       <<" / "<<active->t<<endl
	       <<"  Mom : "<<v_mo<<" : "<<v_mo.Abs2()<<" / "<<mother->t<<endl
	       <<"  Sis : "<<v_si<<" : "<<v_si.Abs2()<<" / "<<sister->t<<endl
	       <<"  nominal e_sis :"<<(1./active->z-1.)*sqrt(sprime/4.)
	       <<" / "<<(1./active->z-1.)*sqrt(sprime/4.)-sister->t/sqrt(4.*sprime)<<endl
	       <<"  Test s' "<<sprime/active->z<<" =?= "
	       <<(v_mo+partner->part->Momentum()).Abs2()<<endl;
    return 0;
  }
  
  mother->part->SetMomentum(v_mo);
  sister->part->SetMomentum(v_si);
  sister->E2 = sqr(v_si[0]);
  if (ResetEnergies(sister)) kink->DoKinematics(sister);
                        else return 0;
  return 1;
};

bool Spacelike_Kinematics::CheckVector(Vec4D vec) {
  if ( (vec.Abs2() > 0) && (vec.Abs2() < 0) ) return 1;
  if (vec[0] < 0) return 1;
  return 0;
};

double Spacelike_Kinematics::BoostInCMS(Tree ** trees,Knot * active, Knot * partner) {
  Vec4D cms = active->part->Momentum()+partner->part->Momentum();

  msg.Debugging()<<"BoostInCMS ... "<<cms<<endl
		 <<"    "<<active->part->Momentum()<<endl
		 <<"  + "<<partner->part->Momentum()<<endl;

  boost     = Poincare(cms);
  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);
  rot       = Poincare(active->part->Momentum(),ZVEC);
  trees[0]->BoRo(rot);
  trees[1]->BoRo(rot);

  return cms.Abs2();
}

Vec4D Spacelike_Kinematics::BoostInLab(Tree ** trees) {
  Knot * init1 = trees[0]->GetInitiator();
  Knot * init2 = trees[1]->GetInitiator();

  double E     = rpa.gen.Ecms();
  Vec4D  cms1  = init1->part->Momentum();
  double E1    = init1->x * E;
  double p1    = sqrt(E1*E1 - cms1.Abs2());

  Vec4D  cms2  = init2->part->Momentum();
  double E2    = init2->x * E;
  double p2    = sqrt(E2*E2 - cms2.Abs2());

  Vec4D  lab1  = Vec4D(E1,0.,0.,p1);
  Vec4D  lab2  = Vec4D(E2,0.,0.,-p2);

  msg.Debugging()<<"BoostInLab ... "<<E1<<", "<<init1->x<<", "<<cms1<<", "<<p1<<endl
		 <<"BoostInLab ... "<<E2<<", "<<init2->x<<", "<<cms2<<", "<<p2<<endl;

  boost        = Poincare(lab1+lab2);

  Vec4D root1  = trees[0]->GetRoot()->part->Momentum();
  Vec4D root2  = trees[1]->GetRoot()->part->Momentum();

  msg.Debugging()<<"Before  :"<<root1<<" : "<<trees[0]->GetRoot()->t<<endl
		 <<"         "<<root2<<" : "<<trees[1]->GetRoot()->t<<endl;

  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);
  rot = Poincare(cms1,lab1);
  trees[0]->BoRo(rot);
  trees[1]->BoRo(rot);

  root1  = trees[0]->GetRoot()->part->Momentum();
  root2  = trees[1]->GetRoot()->part->Momentum();

  msg.Debugging()<<"Finally :"<<root1<<" "<<"         "<<root2<<endl;

  return root1+root2;
}

bool Spacelike_Kinematics::ResetEnergies(Knot * in) {
  msg.Debugging()<<"Spacelike_Kinematics::ResetEnergies(Knot "<<in->kn_no<<" )"<<endl;
  if (in->E2 < in->t) return 0;
  msg.Debugging()<<"    Passed timelike condition : "<<in->E2<<" > "<<in->t<<endl;

  if (in->left) {
    in->left->E2  = in->z*in->z*in->E2;
    in->right->E2 = (1.-in->z)*(1.-in->z)*in->E2;
    if (!ResetEnergies(in->left)) return 0;
    if (!ResetEnergies(in->right)) return 0;
  }
  return 1;
};

bool Spacelike_Kinematics::KinCheck(Knot * active,bool jetveto)
{
  if (!jetveto) return 0;
  if (jf->TwoJets(active->part->Momentum())) return 1;
  return 0;
}

double Spacelike_Kinematics::CalculateMaxT(Knot * active,Knot * partner) {
  if ((!active) || (!partner) || (!active->prev)) {
    msg.Error()<<"Error in Spacelike_Kinematics."<<endl
	       <<"   CalculateMaxT : No knots."<<endl;
    return 0.;
  }
  double t1   = active->t;
  double t2   = partner->t;
  double t3   = active->prev->t;
  double s    = (active->part->Momentum() + partner->part->Momentum()).Abs2();
  double s1   = s           - t2 - t1;
  double s3   = s/active->z - t2 - t3;
  double l1   = s1*s1-4.*t2*t1;
  double l3   = s3*s3-4.*t2*t3;
  if ((l1<0.) || (l3<0.)) return -1.;

  double np1  = sqrt(l1);
  double np3  = sqrt(l3);
  double maxt;
  if (dabs(t2)>rpa.gen.Accu()) maxt = (t1 + t3 + (np1*np3 - s1*s3)/(2.*t2));
                          else maxt = (-(t1/active->z - t3)*(s/(s-t1)-s/(s/active->z-t3))); 
  return maxt;
};






