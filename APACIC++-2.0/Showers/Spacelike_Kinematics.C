#include "Spacelike_Kinematics.H"
#include "Run_Parameter.H"
#include <iomanip>

using std::endl;

using namespace APACIC;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;


void Spacelike_Kinematics::InitKinematics(Tree ** trees,Knot * k1, Knot * k2, int first) 
{
  if ((!k1) || (!k2)) {
    msg.Error()<<"ERROR in Spacelike_Kinematics::InitKinematics : No knots."<<endl;
    return;
  }  

  double t1 = k1->part->Momentum().Abs2();
  double t2 = k2->part->Momentum().Abs2();

  if (first)
  BoostInCMS(trees,k1, k2);

  Vec4D o1 = k1->part->Momentum();
  Vec4D o2 = k2->part->Momentum();

  msg.Tracking() <<"Spacelike_Kinematics::InitKinematics : A "<<endl
   <<"   Vec1 : "<<o1<<" : "<<o1.Abs2()<<" / "<<k1->t<<endl
   <<"   Vec2 : "<<o2<<" : "<<o2.Abs2()<<" / "<<k2->t<<endl
   <<"   S    : "<<(o1+o2).Abs2()<<endl;


  //  rot           = Poincare(Vec4D::ZVEC,o1);
  Vec4D cms     = o1 + o2;
  

  double sprime = cms.Abs2();
  double E1     = (sprime + k1->t - k2->t)/sqrt(4.*sprime);
  double E2     = (sprime - k1->t + k2->t)/sqrt(4.*sprime);
  double pz     = sqrt((sqr(sprime - k1->t - k2->t)-4.*k1->t*k2->t)/(4.*sprime));
  Vec4D v1(E1,0.,0.,pz);
  Vec4D v2(E2,0.,0.,-pz);



  if (first==1) {
    // create boost (for "daughters" of "mother" with unchanged t)
    if (IsEqual(k1 -> t,t1)) {
      double s1 = 4. * sqr(o1[0]);
      double t1 = k1 -> t;
      double sgn1 = E1/dabs(E1);
      double eboo1 = -(dabs(E1) * s1 - pz*sqrt(s1*(s1-4.*t1))); // / -(2. * t1);
      double pboo1 = sgn1 * (pz * s1 - dabs(E1)*sqrt(s1*(s1-4.*t1))); // /-(2. * t1);
      Vec4D b1(eboo1,0.,0.,pboo1);
      boost = Poincare(b1);
      boost.Boost(o1);
      trees[0]->BoRo(boost);
      msg.Tracking() <<"Spacelike_Kinematics::InitKinematics : B "<<endl
		     <<"   Vec1 : "<<o1<<" : "<<o1.Abs2()<<" / "<<k1->t<<endl
		     <<"   Boo1 : "<<b1<<" : "<<b1.Abs2()<<endl;
    }

    if (IsEqual(k2 -> t,t2)) {
      double s2 = 4. * sqr(o2[0]);
      double t2 = k2 -> t;
      double sgn2 = E2/dabs(E2);
      double eboo2 = - (dabs(E2) * s2 - pz*sqrt(s2*(s2-4.*t2))); // /-(2. * t2);
      double pboo2 = - sgn2 * (pz * s2 - dabs(E2)*sqrt(s2*(s2-4.*t2))); // /-/(2. * t2);
      Vec4D b2(eboo2,0.,0.,pboo2);
      boost = Poincare(b2);
      boost.Boost(o2);
      trees[0]->BoRo(boost);
      msg.Tracking() <<"Spacelike_Kinematics::InitKinematics : B "<<endl
		     <<"   Vec2 : "<<o2<<" : "<<o2.Abs2()<<" / "<<k2->t<<endl
		     <<"   Boo2 : "<<b2<<" : "<<b2.Abs2()<<endl;
    }
  }

  k1->part->SetMomentum(v1);
  k2->part->SetMomentum(v2);  

  msg.Tracking() <<"Spacelike_Kinematics::InitKinematics : C"<<endl
		 <<"   Vec1 : "<<v1<<" : "<<v1.Abs2()<<" / "<<k1->t<<endl
		 <<"   Vec2 : "<<v2<<" : "<<v2.Abs2()<<" / "<<k2->t<<endl
		 <<"   S    : "<<(v1+v2).Abs2()<<" / "<<sprime<<endl;
}

bool Spacelike_Kinematics::DoKinematics(Tree ** trees,Knot * active, Knot * partner,int leg, int first) 
{
  if (leg==0) msg.Tracking()<<" Spacelike_Kinematics::DoKinematics( I ";
  else msg.Tracking()<<" Spacelike_Kinematics::DoKinematics( II ";
  msg.Tracking()<<",("<<active->kn_no<<"), <"<<partner->kn_no<<"> );"<<endl;

  if (!active->prev) {
    msg.Tracking()<<"Error Spacelike_Kinematics::DoKinematics : "
	          <<"     No mother for active knot, no kinematics to be constructed"<<endl;
    return 0;
  }

  Vec4D o1 = active->part->Momentum();
  Vec4D o2 = partner->part->Momentum();

  msg.Tracking() <<"Spacelike_Kinematics::DoKinematics : A "<<endl
   <<"   Vec1 : "<<o1<<" : "<<o1.Abs2()<<" / "<<active->t<<endl
   <<"   Vec2 : "<<o2<<" : "<<o2.Abs2()<<" / "<<partner->t<<endl
   <<"   S    : "<<(o1+o2).Abs2()<<endl;


  BoostInCMS(trees,active, partner);

  o1 = active->part->Momentum();
  o2 = partner->part->Momentum();

  msg.Tracking() <<"Spacelike_Kinematics::DoKinematics : B "<<endl
   <<"   Vec1 : "<<o1<<" : "<<o1.Abs2()<<" / "<<active->t<<endl
   <<"   Vec2 : "<<o2<<" : "<<o2.Abs2()<<" / "<<partner->t<<endl
   <<"   S    : "<<(o1+o2).Abs2()<<endl;

  Vec4D cms      = active->part->Momentum()+partner->part->Momentum();
  double sprime  = cms.Abs2();

  Knot * mother  = active->prev;
  Knot * sister  = mother->left;

  double s1      = sprime           - (partner->t) - (active->t);    
  double s3      = sprime/active->z - (partner->t) - (mother->t);
  double np1     = sqrt(s1*s1-4.*(partner->t)*(active->t));
  double np3     = sqrt(s3*s3-4.*(partner->t)*(mother->t));
  double maxt_d2 = CalculateMaxT(active,partner);

  if (maxt_d2 < sister->t) {
    msg.Debugging()<<" if max_d2 < sister->t "<<maxt_d2<<" sister "<<sister->t<<endl;
    return 0;
  }
  
  double E_mo      = 1./(2.*sqrt(sprime)) *
    (sprime/active->z - partner->t + active->t - sister->t);
  double pz_mo     = 1./(2.*active->part->Momentum()[3]) * 
    (s3 - 2. * partner->part->E() * E_mo);
  double pt_mo     = sqrt((maxt_d2-sister->t)/(np1*np1) *
    (0.5*(s1*s3 + np1*np3) + partner->t * (sister->t - active->t - mother->t)));
 

  double test_f4  =((1./active->z-1.)*sprime - sister->t)/(2.*sqrt(sprime));
  msg.Debugging()<<" (spr,z,t4) : "<<sprime<<","<<active->z<<","<<sister->t<<endl;
  msg.Debugging()<<"  e4^2  =="<<sqr(test_f4)<<endl;

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
    msg.Error()<<" (spr,z,t4) : " <<sprime<<","<<active->z<<","<<sister->t<<endl;
    // *AS*    return 0;
  }
  
  mother->part->SetMomentum(v_mo);
  sister->part->SetMomentum(v_si);
  sister->E2 = sqr(v_si[0]);
  if (ResetEnergies(sister)) {
    msg.Debugging()<<" calling finalstate kinematics "<<endl;
    kink->DoKinematics(sister);
  }
  else {
    msg.Debugging()<<" reset energies failed "<<endl;
    return 0;
  }
  return 1;
}

bool Spacelike_Kinematics::CheckVector(Vec4D vec) 
{
  if ( (vec.Abs2() > 0) && (vec.Abs2() < 0) ) return 1;
  if (vec[0] < 0) return 1;
  return 0;
}

double Spacelike_Kinematics::BoostInCMS(Tree ** trees,Knot * active, Knot * partner) 
{
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

Vec4D Spacelike_Kinematics::BoostInLab(Tree ** trees) 
{
  Knot * init1 = trees[0]->GetInitiator();
  Knot * init2 = trees[1]->GetInitiator();

  double E     = rpa.gen.Ecms()/2.;
  Vec4D  cms1  = init1->part->Momentum();
  double E1    = init1->x * E;
  double p1    = sqrt(E1*E1 - cms1.Abs2());

  Vec4D  cms2  = init2->part->Momentum();
  double E2    = init2->x * E;
  double p2    = sqrt(E2*E2 - cms2.Abs2());

  Vec4D  lab1  = Vec4D(E1,0.,0.,p1);
  Vec4D  lab2  = Vec4D(E2,0.,0.,-p2);

  Vec4D  lab   = Vec4D(E1+E2,0.,0.,p2-p1);
  boost        = Poincare(lab);
  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);

  return trees[0]->GetRoot()->part->Momentum() + trees[1]->GetRoot()->part->Momentum();
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
}

bool Spacelike_Kinematics::JetCheck(Knot * active)
{
  Knot * sister=active;
  cout<<" JetCheck ("<<active->kn_no<<")"<<endl;
  cout<<active->part->Momentum()<<endl;
  if (sister->prev) {
    cout<<" JetCheck2 "<<endl;
    sister=sister->prev;
    if (sister->left) {
      cout<<" JetCheck3 "<<endl;

      sister=sister->left;
      if (jf->TwoJets(sister->part->Momentum())) {
	cout<<" JetCheck VETO"<<endl;
  
	msg.Tracking()<<" Spacelike_Kinematics::JetCheck   JETVETO"<<std::endl;
	return 1;
      } 
      cout<<" JetCheck no veto"<<endl;

    }
  }
  return 0;
}

bool Spacelike_Kinematics::KinCheck(Knot * active,bool jetveto)
{
  if (!jetveto) return 0;
  /*
    moved to veto
    if (jf->TwoJets(active->part->Momentum())) {
      msg.Tracking()<<" Spacelike_Kinematics::KinCheck   JETVETO"<<std::endl;
      return 1;
    }
  */
  return 0;
}

double Spacelike_Kinematics::CalculateMaxT(Knot * active,Knot * partner) {
  if ((!active) || (!partner) || (!active->prev)) {
    msg.Error()<<"ERROR: in Spacelike_Kinematics."<<endl
	       <<"   CalculateMaxT : No knots."<<endl;
    return 0.;
  }
  msg.Debugging()<<" in Spacelike_Kinematics::CalculateMaxT( k1: ("<<active->kn_no<<"), k2: ("
      <<partner->kn_no<<"), k3: ("<<active->prev->kn_no<<") )"<<endl;
  double t1   = active->t; 
  double t2   = partner->t;
  double t3   = active->prev->t;
  double s    = (active->part->Momentum() + partner->part->Momentum()).Abs2();
  double s1   = s           - t2 - t1;
  double s3   = s/active->z - t2 - t3;
  double l1   = s1*s1-4.*t2*t1;
  double l3   = s3*s3-4.*t2*t3;
  msg.Debugging()<<"(t1,t2,t3), (s1,s3) : ("<<t1<<","<<t2<<","<<t3<<"), ("<<s1<<","<<s3<<")"<<endl;
  if ((l1<0.) || (l3<0.)) return -1.;

  double np1  = sqrt(l1);
  double np3  = sqrt(l3);
  double maxt;
  if (dabs(t2)>rpa.gen.Accu()) maxt = (t1 + t3 + (np1*np3 - s1*s3)/(2.*t2));
                          else maxt = (-(t1/active->z - t3)*(s/(s-t1)-s/(s/active->z-t3))); 
  return maxt;
}






