#include "Spacelike_Kinematics.H"
#include "Run_Parameter.H"
#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

Spacelike_Kinematics::Spacelike_Kinematics(ATOOLS::Jet_Finder *const jf): 
  p_jf(jf), p_kin(new Timelike_Kinematics(p_jf)) {}

Spacelike_Kinematics::~Spacelike_Kinematics() 
{ 
  if (p_kin) delete p_kin;
}

void Spacelike_Kinematics::InitKinematics(Tree **const trees,Knot *const k1, 
					  Knot *const  k2,const int &first) 
{
  if (k1==NULL || k2==NULL) {
    msg.Error()<<METHOD<<"(..): No knots. Abort."<<std::endl;
    return;
  }  
  if (first) BoostInCMS(trees,k1,k2);
  double t1(k1->part->Momentum().Abs2()), t2(k2->part->Momentum().Abs2());
  Vec4D  o1(k1->part->Momentum()), o2(k2->part->Momentum()), cms(o1+o2);
  double sprime(cms.Abs2()), rtsprime(sqrt(sprime));
  double E1((sprime+k1->t-k2->t)/(2.0*rtsprime)), E2(rtsprime-E1);
  double pz(sqrt(E1*E1-k1->t));
  Vec4D  v1(E1,0.0,0.0,pz), v2(E2,0.0,0.0,-pz);
  static double accu(sqrt(rpa.gen.Accu()));
  Vec4D::SetAccu(accu);
  if (first==1 && (!IsEqual(k1->t,t1)||!IsEqual(k2->t,t2))) {
    if (IsEqual(k1->t,t1)) {
      double sgnt(k1->t<0.0?-1.0:1.0), p1(sqrt(o1[0]*o1[0]-k1->t));
      Vec4D b1(sgnt*(E1*o1[0]-pz*p1),0.0,0.0,-sgnt*(pz*o1[0]-E1*p1));
      m_boost=Poincare(b1);
      m_boost.Boost(o1);
      if (!(o1==v1)) {
	msg.Error()<<METHOD<<"(..): Four momentum not conserved.\n"
		   <<"  p_miss  = "<<(v1-o1)<<"\n"
		   <<"  p_old   = "<<o1<<" "<<o1.Abs2()<<" <- "<<k1->t<<"\n"
		   <<"  p_new   = "<<v1<<" "<<v1.Abs2()<<" <- "<<k1->t<<"\n"
		   <<"  p_boost = "<<b1<<" "<<b1.Abs2()<<std::endl;
      }
      trees[0]->BoRo(m_boost);
    }
    if (IsEqual(k2->t,t2)) {
      double sgnt(k2->t<0.0?-1.0:1.0), p2(sqrt(o2[0]*o2[0]-k2->t));
      Vec4D b2(sgnt*(E2*o2[0]-pz*p2),0.0,0.0,sgnt*(pz*o2[0]-E2*p2));
      m_boost=Poincare(b2);
      m_boost.Boost(o2);
      if (!(o2==v2)) {
	msg.Error()<<METHOD<<"(..): Four momentum not conserved.\n"
		   <<"  p_miss  = "<<(v2-o2)<<"\n"
		   <<"  p_old   = "<<o2<<" "<<o2.Abs2()<<" <- "<<k2->t<<"\n"
		   <<"  p_new   = "<<v2<<" "<<v2.Abs2()<<" <- "<<k2->t<<"\n"
		   <<"  p_boost = "<<b2<<" "<<b2.Abs2()<<std::endl;
      }
      trees[1]->BoRo(m_boost);
    }
  }
  Vec4D::ResetAccu();
  k1->part->SetMomentum(v1);
  k2->part->SetMomentum(v2);  
}

bool Spacelike_Kinematics::DoKinematics(Tree **const trees,Knot *const active,
					Knot *const partner,const int &leg,
					const int &first,const bool test) 
{
  msg_Debugging()<<METHOD<<"("<<active->kn_no<<","<<partner->kn_no
		 <<","<<leg<<"): {\n";
  msg_Indent();
  if (active->prev==NULL) {
    msg.Error()<<METHOD<<"(..): No mother. Abort."<<std::endl;
    msg_Debugging()<<"}\n";
    return false;
  }
  Knot *mother(active->prev), *sister(mother->left);
  int mode(0);
  if (mother->part->Info()=='H') {
    // momenta known (keep phi)   "B"
    mode=1;
    if (mother->prev) if (mother->prev->part->Info()=='H') {
      // boost mother tree appropriate (after momentum determination)  "C1"
      mode+=2;
    }
    if (sister->left) if (sister->left->part->Momentum()[1]!=0.) {
      // if (sister->left->part->Info()=='H') {
      // boost sister tree appropriate (after momentum determination)  "C2"
      mode+=4;
    }
  }
  if (sister->left && mode==0) if (sister->left->part->Momentum()[1]!=0.) {
    mode=5;
  }
  if (test && mode!=5) {
    msg_Debugging()<<"}\n";
    return true;
  }
  msg_Debugging()<<mother->part->Info()<<","
		 <<(mother->prev?mother->prev->part->Info():'-')<<","
		 <<(sister->left?sister->left->part->Momentum()[1]:1.e37)
		 <<" -> mode = "<<mode<<"\n";
  BoostInCMS(trees,active, partner);
  if (mode!=7) {
    Vec4D o1(active->part->Momentum()), o2(partner->part->Momentum()), cms(o1+o2);
    double sprime(cms.Abs2()), s3(sprime/active->z-partner->t-mother->t);
    double maxt_d2(CalculateMaxT(active,partner));
    if (maxt_d2<sister->t) return false;
    double E_mo(1./(2.*sqrt(sprime))*(sprime/active->z-partner->t+active->t-sister->t));
    double pz_mo(1./(2.*active->part->Momentum()[3])*(s3-2.*partner->part->E()*E_mo));
    double pt_mo(sqrt(sqr(E_mo)-sqr(pz_mo)-mother->t)), cph(cos(active->phi)), sph(sin(active->phi));
    if (mode>0) {
      // determine phi
      Vec4D p_mo= mother->part->Momentum();
      double pt=sqrt(sqr(p_mo[1]) + sqr(p_mo[2]));
      if (pt!=0.0) {
	cph=p_mo[2]/pt;
	sph=p_mo[1]/pt;
      }
      else {
	cph=1.0;
	sph=0.0;
      }
      active->phi=acos(cph);
    }
    Vec4D v_mo(E_mo,sph*pt_mo,cph*pt_mo,pz_mo);
    Vec4D v_si(v_mo + (-1.)*active->part->Momentum());
    if (active->part->Momentum().Nan() || partner->part->Momentum().Nan() ||
	v_mo.Nan() || v_si.Nan()) {
      msg.Error()<<METHOD<<"(..): Error. Bad vectors.\n"
		 <<"  Act: "<<active->part->Momentum()<<" "
		 <<active->part->Momentum().Abs2()<<" / "<<active->t<<"\n"
		 <<"  Mom: "<<v_mo<<" "<<v_mo.Abs2()<<" / "<<mother->t<<"\n"
		 <<"  Sis: "<<v_si<<" "<<v_si.Abs2()<<" / "<<sister->t<<"\n"
		 <<"  nom E_sis: "<<(1./active->z-1.)*sqrt(sprime/4.)
		 <<" / "<<((1./active->z-1.)*sqrt(sprime/4.)-
			   sister->t/sqrt(4.*sprime))<<"\n"
		 <<"  test s' "<<sprime/active->z<<" =?= "
		 <<(v_mo+partner->part->Momentum()).Abs2()<<"\n"
		 <<" (spr,z,t4): " <<sprime<<", "<<active->z<<", "
		 <<sister->t<<std::endl;
    }
    BoostPartial(mode,mother,sister,v_mo,v_si);
  }
  if (test && mode==5) {
    BoostFromCMS(trees);
    msg_Debugging()<<"}\n";
    return true;
  }
  if (ResetEnergies(sister)) p_kin->DoKinematics(sister);
  else {
    msg_Debugging()<<"}\n";
    return false;
  }
  msg_Debugging()<<"}\n";
  return true;
}


void Spacelike_Kinematics::BoostPartial(const int mode,Knot *const si,
					const Vec4D &v_si) 
{
  msg_Debugging()<<METHOD<<"("<<mode<<","<<si->kn_no<<","<<v_si<<"): {\n";
  msg_Indent();
  Vec4D p_si(si->part->Momentum());
  double E11(p_si[0]), t1(p_si.Abs2());
  static double accu(sqrt(rpa.gen.Accu()));
  if (dabs((t1-v_si.Abs2())/t1)>accu) {
    msg.Error()<<METHOD<<"(..): Mass deviation. t1 = "<<t1
	       <<" t1-t1' = "<<t1-v_si.Abs2()<<std::endl;
  }
  // rotation
  Poincare rot(p_si,v_si);
  rot.Rotate(p_si);
  Vec3D  np_si(p_si/p_si.PSpat()), nv_si(v_si/v_si.PSpat());
  double E1(v_si[0]), pz(v_si.PSpat());
  double sgnt(t1<0.0?-1.0:1.0);
  // boost
  Vec4D b1(sgnt*(E1*E11-pz*sqrt(E11*E11-t1)),
	   -sgnt*(pz*E11-E1*sqrt(E11*E11-t1))*np_si);
  m_boost=Poincare(b1);
  m_boost.Boost(p_si);
  // apply rotation and boost on tree rest (sequence depending on mode)
  msg_Debugging()<<"mode = "<<mode<<"\n";
  if (mode==3) RoBoIni(si,rot,m_boost);
  else if (mode==5) RoBoFin(si,rot,m_boost);
  else msg.Error()<<METHOD<<"(..): Error."<<std::endl;
  msg_Debugging()<<"}\n";
}


void Spacelike_Kinematics::RoBoIni(Knot *const k,Poincare &rot,Poincare &boost) 
{
  if (k==NULL) return;
  msg_Debugging()<<METHOD<<"("<<k->kn_no<<"):\n";
  Vec4D  p(k->part->Momentum());
  rot.Rotate(p);
  boost.Boost(p);
  k->part->SetMomentum(p);
  if (k->prev) {
    RoBoIni(k->prev,rot,boost);
    RoBoFin(k->prev->left,rot,boost);
  }
}

void Spacelike_Kinematics::RoBoFin(Knot *const k,Poincare &rot,Poincare &boost) 
{
  if (k==NULL) return;
  msg_Debugging()<<METHOD<<"("<<k->kn_no<<"):\n";
  Vec4D  p(k->part->Momentum());
  rot.Rotate(p);
  boost.Boost(p);
  k->part->SetMomentum(p);
  RoBoFin(k->left,rot,boost);
  RoBoFin(k->right,rot,boost);
}

void Spacelike_Kinematics::BoostPartial(const int mode,
					Knot *const mo,Knot *const si,
					const Vec4D &v_mo,const Vec4D &v_si) 
{
  msg_Debugging()<<METHOD<<"("<<mode<<","<<mo->kn_no<<","
		 <<si->kn_no<<","<<v_mo<<","<<v_si<<"): {\n";
  msg_Indent();
  msg_Debugging()<<"mode = "<<mode<<"\n";
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
    if (si->t!=si->tout) {
      double tnew(v_si.Abs2());
      static double accu(sqrt(rpa.gen.Accu()));
      if (dabs(tnew/si->t-1.0)>accu)  
	msg.Error()<<METHOD<<"(..): Relative mass deviation "
		   <<(tnew/si->t-1.0)<<". Reset t."<<std::endl;
      si->t=tnew;
    }
    mo->part->SetMomentum(v_mo);
    si->part->SetMomentum(v_si);
    si->E2 = sqr(v_si[0]);
  }
  msg_Debugging()<<"}\n";
}

double Spacelike_Kinematics::BoostInCMS(Tree **const trees,
					Knot *const k1,Knot *const k2) 
{
  Vec4D cms(k1->part->Momentum()+k2->part->Momentum());
  m_boost=Poincare(cms);
  trees[0]->BoRo(m_boost);
  trees[1]->BoRo(m_boost);
  m_rot=Poincare(k1->part->Momentum(),Vec4D::ZVEC);
  trees[0]->BoRo(m_rot);
  trees[1]->BoRo(m_rot);
  return cms.Abs2();
}

double Spacelike_Kinematics::BoostFromCMS(Tree **const trees) 
{
  m_rot.Invert();
  trees[0]->BoRo(m_rot);
  trees[1]->BoRo(m_rot);
  m_boost.Invert();
  trees[0]->BoRo(m_boost);
  trees[1]->BoRo(m_boost);
  m_rot.Invert();
  m_boost.Invert();
  return 0.;
}

Vec4D Spacelike_Kinematics::BoostInLab(Tree **const trees) 
{
  double x1(trees[0]->GetInitiator()->x), x2(trees[1]->GetInitiator()->x);
  // Only for massless initiators.
  Vec4D  lab(Vec4D(x1+x2,0.,0.,x2-x1));
  m_boost = Poincare(lab);
  trees[0]->BoRo(m_boost);
  trees[1]->BoRo(m_boost);
  return trees[0]->GetRoot()->part->Momentum()+ 
    trees[1]->GetRoot()->part->Momentum();
}

bool Spacelike_Kinematics::ResetEnergies(Knot *const in) 
{
  if (in->E2 < in->t) return 0;
  if (in->left) {
    if (in->part->Info()=='H' && in->left->part->Info()=='H') {
      // update z:
      in->z=in->left->part->Momentum()[0]/in->part->Momentum()[0];
    }
    in->left->E2  = in->z*in->z*in->E2;
    in->right->E2 = (1.-in->z)*(1.-in->z)*in->E2;
    if (!ResetEnergies(in->left)) return 0;
    if (!ResetEnergies(in->right)) return 0;
  }
  return 1;
}

double Spacelike_Kinematics::CalculateMaxT(Knot *const active,
					   Knot *const partner) 
{
  if (active==NULL || partner==NULL || active->prev==NULL) {
    msg.Error()<<"ERROR: in Spacelike_Kinematics."<<std::endl
	       <<"   CalculateMaxT : No knots."<<std::endl;
    return 0.;
  }
  double t1(active->t), t2(partner->t), t3(active->prev->t);
  double s((active->part->Momentum() + partner->part->Momentum()).Abs2());
  double s1(s-t2-t1), s3(s/active->z-t2-t3);
  double l1(s1*s1-4.*t2*t1), l3(s3*s3-4.*t2*t3);
  if (l1<0. || l3<0.) return -1.;
  double np1(sqrt(l1)), np3(sqrt(l3));
  double maxt0(-(t1/active->z - t3)*(s/(s-t1)-s/(s/active->z-t3))); 
  double qa(- s1*s3), qb(np1*np3), q(qa + qb);
  if (dabs(q/qb)<1.e-8)  {
    double z(active->z), sprime(s), sprime2(sqr(sprime)), sprime3(sprime*sprime2), sprime4(sqr(sprime2));
    double t12(sqr(t1)), t13(t1*t12), t32(sqr(t3)), t33(t3*t32), z2(sqr(z)), z3(z*z2);
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

void Spacelike_Kinematics::ResetMomenta(Knot *const k, Tree *const tree, 
					const int &mode) 
{
  if (mode==0) {
    Knot *mo(k->prev), *si(mo->left);
    // kill not-ME daughters of sister
    ResetMomenta(si,tree,1);
    // kill not-ME grand mother
    ResetMomenta(mo,tree,2);
    k->stat=1;
  }
  else if (mode==1) {
    Knot *si(k);
    if (si->left) {
      if (si->left->part->Info()!='H') {
	si->left=NULL;
	si->right=NULL;
	si->stat=3;
      }
      else {
	ResetMomenta(si->left,tree,1);
	ResetMomenta(si->right,tree,1);
      }
    }
    else {
      si->stat=3;
      si->t=0.;
      si->tmax=0.;
    }
  }
  else if (mode==2) {
    Knot *mo=k;
    if (mo->prev) {
      if (mo->prev->part->Info()!='H') {
	mo->prev=NULL;
	mo->stat=3;
      }
    }
    else {
      mo->stat=3;
    }
  }
}
