#include "Spacelike_Kinematics.H"
#include "Run_Parameter.H"
#include <iomanip>

using std::endl;

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

Spacelike_Kinematics::Spacelike_Kinematics(double pt2minFS, Data_Read * const dataread) {
  double ycut   = ATOOLS::rpa.gen.Ycut();
  m_type = 1;
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton()) {
    m_type = 1;
  }
  else if ((!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())) {
    m_type = 4;
  }
  else {
    msg.Error()<<"ERROR in Spacelike_Kinematics : "<<std::endl
	       <<"   DIS is not yet implemented in the jetfinder, continue with hh-mode."<<endl;
    m_type = 4;
  }
  jf   = new ATOOLS::Jet_Finder(ycut,m_type);
  jf->SetDeltaR(rpa.gen.DeltaR()); 
  kink = new Timelike_Kinematics(pt2minFS,dataread);
}


void Spacelike_Kinematics::InitKinematics(Tree ** trees,Knot * k1, Knot * k2, int first) 
{
  if ((!k1) || (!k2)) {
    msg.Error()<<"ERROR in Spacelike_Kinematics::InitKinematics : No knots."<<endl;
    return;
  }  

  
  double t1 = k1->part->Momentum().Abs2();
  double t2 = k2->part->Momentum().Abs2();

  if (first) BoostInCMS(trees,k1, k2);

  Vec4D o1 = k1->part->Momentum();
  Vec4D o2 = k2->part->Momentum();
  Vec4D cms     = o1 + o2;
  

  double sprime = cms.Abs2();
  double E1     = (sprime + k1->t - k2->t)/sqrt(4.*sprime);
  double E2     = (sprime - k1->t + k2->t)/sqrt(4.*sprime);
  double pz     = sqrt((sqr(sprime - k1->t - k2->t)-4.*k1->t*k2->t)/(4.*sprime));
  Vec4D v1(E1,0.,0.,pz);
  Vec4D v2(E2,0.,0.,-pz);


  bool error =0 ;
  if (first==1 && (!IsEqual(k1 -> t,t1)||!IsEqual(k2 -> t,t2))) {
    // create boost (for "daughters" of "mother" with unchanged t)
    if (IsEqual(k1 -> t,t1)) {
      double s1 = 4. * sqr(o1[0]);
      double t1 = k1 -> t;
      double sgn2 = 1.;
      if (t1<0.) sgn2 = -1.;
      double sgn1 = E1/dabs(E1);
      double eboo1 = sgn2*(dabs(E1) * s1 - pz*sqrt(s1*(s1-4.*t1))); // / -(2. * t1);
      double pboo1 = -sgn2*sgn1 * (pz * s1 - dabs(E1)*sqrt(s1*(s1-4.*t1))); // /-(2. * t1);
      Vec4D b1(eboo1,0.,0.,pboo1);
      boost = Poincare(b1);
      boost.Boost(o1);
      if (!(o1==v1)) {
	error=1;
	msg.Out()<<"WARNING in  Spacelike_Kinematics::InitKinematics : Mismatch."<<endl;
	msg.Out().precision(12);
      }
      trees[0]->BoRo(boost);
      if (error) {
	msg.Out() <<"Spacelike_Kinematics::InitKinematics : B "<<endl
		  <<"   Vec1 : "<<o1<<" : "<<o1.Abs2()<<" / "<<k1->t<<endl
		  <<"   Boo1 : "<<b1<<" : "<<b1.Abs2()<<endl;
      }
    }

    if (IsEqual(k2 -> t,t2)) {
      double s2 = 4. * sqr(o2[0]);
      double t2 = k2 -> t;

      double sgn2 = 1.;
      if (t2<0.) sgn2 = -1.;
      double eboo2 = sgn2 * (dabs(E2) * s2 - pz*sqrt(s2*(s2-4.*t2))); // /-(2. * t2);
      double pboo2 = sgn2 * (pz * s2 - dabs(E2)*sqrt(s2*(s2-4.*t2))); // /-/(2. * t2);
      Vec4D b2(eboo2,0.,0.,pboo2);
      boost = Poincare(b2);
      boost.Boost(o2);
      trees[1]->BoRo(boost);
      if (!(o2==v2)) {
	error=1;
	msg.Out()<<"WARNING in Spacelike_Kinematics::InitKinematics : Mismatch."<<endl;
	msg.Out().precision(12);
      }
      if (error) {
	msg.Out() <<"Spacelike_Kinematics::InitKinematics : B "<<endl
		  <<"   Vec2 : "<<o2<<" : "<<o2.Abs2()<<" / "<<k2->t<<endl
		  <<"   Boo2 : "<<b2<<" : "<<b2.Abs2()<<endl;
      }
    }
  }

  k1->part->SetMomentum(v1);
  k2->part->SetMomentum(v2);  

  if (error) {
    msg.Out()<<"Spacelike_Kinematics::InitKinematics : C"<<endl
	     <<"   Vec1 : "<<v1<<" : "<<v1.Abs2()<<" / "<<k1->t<<endl
	     <<"   Vec2 : "<<v2<<" : "<<v2.Abs2()<<" / "<<k2->t<<endl
	     <<"   S    : "<<(v1+v2).Abs2()<<" / "<<sprime<<endl;
  }


}

bool Spacelike_Kinematics::DoKinematics(Tree ** trees,Knot * active, Knot * partner,int leg, int first, bool test) 
{
  if (!active->prev) {
    msg.Tracking()<<"Error Spacelike_Kinematics::DoKinematics : "
	          <<"     No mother for active knot, no kinematics to be constructed"<<endl;
    return 0;
  }


  Knot * mother  = active->prev;
  Knot * sister  = mother->left;


  int mode=0;
  if (mother->part->Info()=='H') {
    // momenta known (keep phi)   "B"
    mode=1;
    if (mother->prev) if (mother->prev->part->Info()=='H') {
      // boost mother tree appropriate (after momentum determination)  "C1"
      mode+=2;
    }
    if (sister->left) if (sister->left->part->Info()=='H') {
      // boost sister tree appropriate (after momentum determination)  "C2"
      mode+=4;
    }
  }


  if (test && mode!=5) return 1;

  BoostInCMS(trees,active, partner);


  if (mode!=7) {

    Vec4D o1 = active->part->Momentum();
    Vec4D o2 = partner->part->Momentum();

    Vec4D cms      = active->part->Momentum()+partner->part->Momentum();
    double sprime  = cms.Abs2();


    double s1      = sprime           - (partner->t) - (active->t);    
    double s3      = sprime/active->z - (partner->t) - (mother->t);
    double np1     = sqrt(s1*s1-4.*(partner->t)*(active->t));
    double np3     = sqrt(s3*s3-4.*(partner->t)*(mother->t));
    double maxt_d2 = CalculateMaxT(active,partner);

    if (maxt_d2 < sister->t) {
      return 0;
    }
  
    double E_mo      = 1./(2.*sqrt(sprime)) *
      (sprime/active->z - partner->t + active->t - sister->t);
    double pz_mo     = 1./(2.*active->part->Momentum()[3]) * 
      (s3 - 2. * partner->part->E() * E_mo);
    double pt_mo     = sqrt((maxt_d2-sister->t)/(np1*np1) *
			    (0.5*(s1*s3 + np1*np3) + partner->t * (sister->t - active->t - mother->t)));
 
    double cph=cos(active->phi), sph=sin(active->phi);

    if (mode>0) {
      // determine phi
      Vec4D p_mo= mother->part->Momentum();
      double pt=sqrt(sqr(p_mo[1]) + sqr(p_mo[2]));
      cph=p_mo[2]/pt;
      sph=p_mo[1]/pt;
      active->phi=acos(cph);
    }

  
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
		 <<(v_mo+partner->part->Momentum()).Abs2()<<endl
		 <<" (spr,z,t4) : " <<sprime<<","<<active->z<<","<<sister->t<<endl;
    }
  
    BoostPartial(mode,mother,sister,v_mo,v_si);
  }

  if (test && mode==5) {
    BoostFromCMS(trees);
    return 1;
  }


  if (ResetEnergies(sister)) {
    kink->DoKinematics(sister);
  }
  else {
    return 0;
  }
  return 1;
}


void Spacelike_Kinematics::BoostPartial(const int mode, Knot * si, const Vec4D & v_si) {
  Vec4D p_si=si->part->Momentum();
  double s1 = 4. * sqr(p_si[0]);
  double t1 = p_si.Abs2();

  if (dabs((t1-v_si.Abs2())/t1)>1.e-7) {
    msg.Error()<<" Warning: mass deviation in Spacelike_Kinematics::BoostPartial "<<std::endl
	       <<" t1="<<t1<<"   t1'="<<v_si.Abs2()<<"  ("<<t1-v_si.Abs2()<<std::endl;
  }

  // rotation
  Poincare rot(p_si,v_si);
  rot.Rotate(p_si);

  Vec3D  np_si(p_si);
  double lp_si=np_si.Abs();
  np_si=np_si*(1./lp_si);

  Vec3D  nv_si(v_si);
  double lv_si=nv_si.Abs();
  nv_si=nv_si*(1./lv_si);

  double E1     = v_si[0];
  double pz     = lv_si;

  double sgn2 = 1.;
  if (t1<0.) sgn2 = -1.;
  double sgn1 = E1/dabs(E1);
  double eboo1 = sgn2*(dabs(E1) * s1 - pz*sqrt(s1*(s1-4.*t1))); // / -(2. * t1);
  double pboo1 = -sgn2*sgn1 * (pz * s1 - dabs(E1)*sqrt(s1*(s1-4.*t1))); // /-(2. * t1);

  // boost
  Vec4D b1(eboo1,pboo1*np_si);
  boost = Poincare(b1);
  boost.Boost(p_si);


  // apply rotation and boost on tree rest (sequence depending on mode)
  if (mode==3) RoBoIni(si,rot,boost);
  else if (mode==5) RoBoFin(si,rot,boost);
  else {
    msg.Error()<<" Error in Spacelike_Kinematics::BoostPartial "<<endl;
  }
}


void Spacelike_Kinematics::RoBoIni(Knot * k, Poincare & rot, Poincare & boost) 
{
  if (!k) return;
  Vec4D  p = k->part->Momentum();
  rot.Rotate(p);
  boost.Boost(p);
  k->part->SetMomentum(p);
  if (k->prev) {
    RoBoIni(k->prev,rot,boost);
    RoBoFin(k->prev->left,rot,boost);
  }
}

void Spacelike_Kinematics::RoBoFin(Knot * k, Poincare & rot, Poincare & boost) 
{
  if (!k) return;
  Vec4D  p = k->part->Momentum();
  rot.Rotate(p);
  boost.Boost(p);
  k->part->SetMomentum(p);
  RoBoFin(k->left,rot,boost);
  RoBoFin(k->right,rot,boost);
}



void Spacelike_Kinematics::BoostPartial(const int mode, Knot * mo, Knot * si, const Vec4D & v_mo, const Vec4D & v_si) 
{
  if (mode==3) {
    // boost mother and the rest
    BoostPartial(mode,mo,v_mo);
    si->part->SetMomentum(v_si);
    si->E2 = sqr(v_si[0]);
  } 
  else if (mode==5) {
    // boost sister and the rest
    BoostPartial(mode,si,v_si);
      
    mo->part->SetMomentum(v_mo);
    si->E2 = sqr(si->part->Momentum()[0]);
  }
  else {
    mo->part->SetMomentum(v_mo);
    si->part->SetMomentum(v_si);
    si->E2 = sqr(v_si[0]);
  }
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

  boost     = Poincare(cms);
  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);
  rot       = Poincare(active->part->Momentum(),Vec4D::ZVEC);
  trees[0]->BoRo(rot);
  trees[1]->BoRo(rot);
  return cms.Abs2();
}

double Spacelike_Kinematics::BoostFromCMS(Tree ** trees) 
{
  rot.Invert();
  trees[0]->BoRo(rot);
  trees[1]->BoRo(rot);
  boost.Invert();
  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);
  rot.Invert();
  boost.Invert();
  return 0.;
}

Vec4D Spacelike_Kinematics::BoostInLab(Tree ** trees) 
{
  double x1 = trees[0]->GetInitiator()->x;
  double x2 = trees[1]->GetInitiator()->x;
  // Only for massless initiators.
  Vec4D  lab   = Vec4D(x1+x2,0.,0.,x2-x1);
  boost        = Poincare(lab);
  trees[0]->BoRo(boost);
  trees[1]->BoRo(boost);

  return trees[0]->GetRoot()->part->Momentum() + trees[1]->GetRoot()->part->Momentum();
}

bool Spacelike_Kinematics::ResetEnergies(Knot * in) {
  if (in->E2 < in->t) return 0;

  if (in->left) {
//     if (in->part->Info()=='H' && in->left->part->Info()=='H') {
//       // update z:
//       in->z=in->left->part->Momentum()[0]/in->part->Momentum()[0];
//     }

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
  if (sister->prev) {
    sister=sister->prev;
    if (sister->left) {
      sister=sister->left;
      if (jf->TwoJets(sister->part->Momentum())) return 1;
    }
  }
  return 0;
}


bool Spacelike_Kinematics::JetVeto(const Vec4D & v)
{
  if (jf->TwoJets(v)) return 1;
  return 0;
}


bool Spacelike_Kinematics::JetVeto(Knot * k1, Knot * k2) 
{
  Vec4D p4=k1->prev->left->part->Momentum();
  Vec4D p1=k1->part->Momentum();
  Vec4D p3=k1->prev->part->Momentum();
  Vec4D p2=k2->part->Momentum();
  double pt2=sqr(p4[1])+sqr(p4[2]);

  Vec4D cms = p3+p2;
  Poincare boost(cms);
  boost.Boost(p1);
  boost.Boost(p2);
  boost.Boost(p3);
  boost.Boost(p4);
  Poincare rot(p3,Vec4D::ZVEC);
  rot.Rotate(p1);
  rot.Rotate(p2);
  rot.Rotate(p3);
  rot.Rotate(p4);
  pt2=sqr(p4[1])+sqr(p4[2]);

  if (jf->TwoJets(p4)) return 1;
  return 0;
}



bool Spacelike_Kinematics::KinCheck(Knot * active,bool jetveto)
{
  if (!jetveto) return 0;
  return 0;
}

double Spacelike_Kinematics::CalculateMaxT(Knot * active,Knot * partner) {
  if ((!active) || (!partner) || (!active->prev)) {
    msg.Error()<<"ERROR: in Spacelike_Kinematics."<<endl
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

  double maxt0 =(-(t1/active->z - t3)*(s/(s-t1)-s/(s/active->z-t3))); 

  double qa =  - s1*s3;
  double qb = np1*np3;
  double q = qa + qb;
  if (dabs(q/qb)<1.e-8)  {
    double z = active->z;
    double sprime  = s;
    double sprime2 = sqr(sprime);
    double sprime3 = sprime*sprime2;
    double sprime4 = sqr(sprime2);
    double t12 = sqr(t1);
    double t13 = t1*t12;
    double t32 = sqr(t3);
    double t33 = t3*t32;
    double z2 = sqr(z);
    double z3 = z*z2;

    double maxt1 = (sprime*t2*(-(t1*(t1 - t3)*t3*z3*(t12 - t32*z)) + sprime4*(-1 + z)*(t1 - t3*z2) + 
		2*sprime*t1*t3*z2*(-2*t1*t3*z - t32*(-2 + z)*z + t12*(-1 + 2*z)) + 
		sprime2*z*(t13 + t12*t3*(5 - 6*z)*z + t33*z3 + t1*t32*z*(-6 + 5*z)) - 
		2*sprime3*z*(t12 + t32*z2 - 2*t1*t3*(1 - z + z2))))/
      (pow((sprime - t1)*(sprime - t3*z),3.)*z);
    return maxt0+maxt1;
  }

  if (dabs(t2)>rpa.gen.Accu()) return  (t1 + t3 + (np1*np3 - s1*s3)/(2.*t2));
  return maxt0;
}

void Spacelike_Kinematics::ResetMomenta(Knot* k, Tree* tree, int mode) 
{
  if (mode==0) {
    Knot *mo=k->prev;
    Knot *si=mo->left;

    // kill not-ME daughters of sister
    ResetMomenta(si,tree,1);
    // kill not-ME grand mother
    ResetMomenta(mo,tree,2);
    k->stat=1;
  }
  else if (mode==1) {
    Knot *si=k;
    if (si->left) {
      if (si->left->part->Info()!='H') {
	si->left=0;
	si->right=0;
	si->stat=3;
      }
      else {
	ResetMomenta(si->left,tree,1);
	ResetMomenta(si->right,tree,1);
      }
    }
    else {
      si->stat=3;
      si->t=0;
      si->tmax=0;
    }
  }
  else if (mode==2) {
    Knot *mo=k;
    if (mo->prev) {
      if (mo->prev->part->Info()!='H') {
	mo->prev=0;
	mo->stat=3;
      }
    }
    else {
      mo->stat=3;
    }
  }

}
