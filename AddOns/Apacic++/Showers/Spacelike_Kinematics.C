#include "AddOns/Apacic++/Showers/Spacelike_Kinematics.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;
using namespace std;

Spacelike_Kinematics::Spacelike_Kinematics(): 
  p_kin(new Timelike_Kinematics()) {}

Spacelike_Kinematics::~Spacelike_Kinematics() 
{ 
  if (p_kin) delete p_kin;
}

int Spacelike_Kinematics::
InitKinematics(Tree **const trees,Tree *const tree,const int tree1,
	       Knot *const k1,Knot *const  k2) 
{
  msg_Debugging()<<METHOD<<"("<<k1->kn_no<<","<<k2->kn_no<<") {\n";
  msg_Indent();
  if (k1==NULL || k2==NULL) {
    msg_Error()<<METHOD<<"(..): No knots. Abort."<<std::endl;
    return -1;
  }
  if (!k1->shower) {
    double s12((k1->part->Momentum()+k2->part->Momentum()).Abs2());
    double s32((k1->prev->part->Momentum()+k2->part->Momentum()).Abs2());
    k1->z=s12/s32;
  }
  if (!k2->shower) {
    double s21((k2->part->Momentum()+k1->part->Momentum()).Abs2());
    double s51((k2->prev->part->Momentum()+k1->part->Momentum()).Abs2());
    k2->z=s21/s51;
  }
  if (k1->shower && k1->t<k1->tout) {
    double maxt(CalculateMaxT(k1,k2,false));
    double tsi(k1->prev->left->stat>0?
	       k1->prev->left->tout:k1->prev->left->t);
    if (maxt<tsi) return 0;
  }
  if (k2->shower && k2->t<k2->tout) {
    double maxt(CalculateMaxT(k2,k1,false));
    double tsi(k2->prev->left->stat>0?
	       k2->prev->left->tout:k2->prev->left->t);
    if (maxt<tsi) return 0;
  }
  double dir(BoostInCMS(trees,tree,tree1,k1,k2));
  Vec4D  o1(k1->part->Momentum()), o2(k2->part->Momentum()), cms(o1+o2);
  double t1(o1.Abs2()), t2(o2.Abs2());
  double sprime(cms.Abs2()), rtsprime(sqrt(sprime));
  double E1((sprime+k1->t-k2->t)/(2.0*rtsprime)), E2(rtsprime-E1);
  double pz(sqrt(E1*E1-k1->t));
  Vec4D  v1(E1,0.0,0.0,dir*pz), v2(E2,0.0,0.0,-dir*pz);
  msg_Debugging()<<"t_{"<<k1->kn_no<<"} = "<<k1->t<<" / "<<t1
		 <<" ("<<(k1->t/t1-1.0)<<") "<<k1->prev<<", t_{"
		 <<k2->kn_no<<"} = "<<k2->t<<" / "<<t2
		 <<" ("<<(k2->t/t2-1.0)<<") "<<k2->prev<<"\n";
  if (dabs(k1->t/t1-1.0)<rpa.gen.Accu() && k1->prev!=NULL) {
    msg_Debugging()<<"boost attached "<<k1->kn_no<<"\n";
    double sgnt(k1->t<0.0?-1.0:1.0), p1(sqrt(o1[0]*o1[0]-k1->t));
    Vec4D b1(sgnt*(E1*o1[0]-pz*p1),0.0,0.0,dir*sgnt*(E1*p1-pz*o1[0]));
    m_boost=Poincare(b1);
    m_boost.Boost(o1);
    if (!IsEqual(o1,v1,sqrt(rpa.gen.Accu()))) {
      msg_Error()<<METHOD<<"(..): Four momentum not conserved on tree 1.\n"
		 <<"  p_miss  = "<<(v1-o1)<<"\n"
		 <<"  p_old   = "<<o1<<" "<<o1.Abs2()<<" <- "<<k1->t<<"\n"
		 <<"  p_new   = "<<v1<<" "<<v1.Abs2()<<" <- "<<k1->t<<"\n"
		 <<"  p_boost = "<<b1<<" "<<b1.Abs2()<<std::endl;
      return -1;
    }
    trees[tree1]->BoRo(m_boost);
  }
  if (dabs(k2->t/t2-1.0)<rpa.gen.Accu() && k2->prev!=NULL) {
    msg_Debugging()<<"boost attached "<<k2->kn_no<<"\n";
    double sgnt(k2->t<0.0?-1.0:1.0), p2(sqrt(o2[0]*o2[0]-k2->t));
    Vec4D b2(sgnt*(E2*o2[0]-pz*p2),0.0,0.0,-dir*sgnt*(E2*p2-pz*o2[0]));
    m_boost=Poincare(b2);
    m_boost.Boost(o2);
    if (!IsEqual(o2,v2,sqrt(rpa.gen.Accu()))) {
      msg_Error()<<METHOD<<"(..): Four momentum not conserved on tree 2.\n"
		 <<"  p_miss  = "<<(v2-o2)<<"\n"
		 <<"  p_old   = "<<o2<<" "<<o2.Abs2()<<" <- "<<k2->t<<"\n"
		 <<"  p_new   = "<<v2<<" "<<v2.Abs2()<<" <- "<<k2->t<<"\n"
		 <<"  p_boost = "<<b2<<" "<<b2.Abs2()<<std::endl;
      return -1;
    }
    trees[1-tree1]->BoRo(m_boost);
  }
  k1->part->SetMomentum(v1);
  k2->part->SetMomentum(v2);  
  msg_Debugging()<<"}\n";
  return 1;
}

bool Spacelike_Kinematics::DoKinematics(Tree **const trees,Tree *const tree,
					const int &leg,
					Knot *const active,Knot *const partner,
					const bool test) 
{
  msg_Debugging()<<METHOD<<"("<<active->kn_no<<","<<partner->kn_no
		 <<","<<leg<<"): {\n";
  msg_Indent();
  if (active->prev==NULL) {
    msg_Error()<<METHOD<<"(..): No mother. Abort."<<std::endl;
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
  BoostInCMS(trees,tree,leg,active,partner);
  if (mode!=7) {
    Vec4D o1(active->part->Momentum()), o2(partner->part->Momentum()), cms(o1+o2);
    double sprime(cms.Abs2()), s3(sprime/active->z-partner->t-mother->t);
    if (s3<0.0) return false;
    double maxt_d2(CalculateMaxT(active,partner));
    double t_si(sister->shower<=1?sister->t:sister->part->Momentum().Abs2());
    if (maxt_d2<t_si) return false;
    double E_mo(1./(2.*sqrt(sprime))*(sprime/active->z-partner->t+active->t-t_si));
    double pz_mo(1./(2.*active->part->Momentum()[3])*(s3-2.*partner->part->E()*E_mo));
    double pt_mo(sqr(E_mo)-sqr(pz_mo)-mother->t), cph(cos(active->phi)), sph(sin(active->phi));
    if (pt_mo<0.0) return false;
    pt_mo=sqrt(pt_mo);
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
      msg_Error()<<METHOD<<"(..): Error. Bad vectors.\n"
		 <<"  Act: "<<active->part->Momentum()<<" "
		 <<active->part->Momentum().Abs2()<<" / "<<active->t<<"\n"
		 <<"  Mom: "<<v_mo<<" "<<v_mo.Abs2()<<" / "<<mother->t<<"\n"
		 <<"  Sis: "<<v_si<<" "<<v_si.Abs2()<<" / "<<t_si
		 <<" ("<<sister->t<<")\n"
		 <<"  nom E_sis: "<<(1./active->z-1.)*sqrt(sprime/4.)
		 <<" / "<<((1./active->z-1.)*sqrt(sprime/4.)-
			   t_si/sqrt(4.*sprime))<<"\n"
		 <<"  test s' "<<sprime/active->z<<" =?= "
		 <<(v_mo+partner->part->Momentum()).Abs2()<<"\n"
		 <<"  (spr,z,t4): " <<sprime<<", "<<active->z<<", "<<t_si<<"\n"
		 <<"  (E3,pz3,pt3,\\sqrt{mo->t},m3): "<<E_mo<<", "<<pz_mo<<", "
		 <<pt_mo<<" "<<sqrt(dabs(mother->t))<<" "<<sqrt(s3)<<std::endl;
    }
    BoostPartial(mode,mother,sister,v_mo,v_si);
  }
  if (test && mode==5) {
    BoostFromCMS(trees,tree);
    msg_Debugging()<<"}\n";
    return true;
  }
  if (ResetEnergies(sister)) {
    static double accu(sqrt(sqrt(Accu())));
    if (sister->part->Info()=='H' && sister->left!=NULL && 
	sister->left->part->Info()=='H' &&
	!IsEqual(sister->t,sister->part->Momentum().Abs2(),accu)) {
      msg_Error()<<METHOD<<"(): Large sister mass deviation: "
		 <<sqrt(sister->part->Momentum().Abs2())<<" vs. "
		 <<sqrt(sister->t)<<" -> "<<
	sister->part->Momentum().Abs2()/sister->t-1.0<<std::endl;
      msg_Debugging()<<"}\n";
      return false;
    }
    p_kin->DoKinematics(sister);
  }
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
    msg_Error()<<METHOD<<"(..): Mass deviation. t1 = "<<t1
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
  if (mode==3) RoBoIni(si,rot,m_boost);
  else if (mode==5) RoBoFin(si,rot,m_boost);
  else msg_Error()<<METHOD<<"(..): Error."<<std::endl;
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
    Vec4D pm(k->prev->part->Momentum());
    Vec4D pd(k->part->Momentum()+k->prev->left->part->Momentum());
    static double accu(sqrt(rpa.gen.Accu()));
    if (!IsEqual(pm,pd,accu)) {
      msg_Error()<<METHOD<<"(..): Four momentum not conserved.\n"
		 <<"  p_miss  = "<<(pm-pd)<<"\n"
		 <<"  p_old   = "<<pm<<" "<<pm.Abs2()<<"\n"
		 <<"  p_new   = "<<pd<<" "<<pd.Abs2()<<std::endl;
    }
  }
}

void Spacelike_Kinematics::RoBoFin(Knot *const k,Poincare &rot,Poincare &boost) 
{
  if (k==NULL) return;
  Vec4D  p(k->part->Momentum());
  rot.Rotate(p);
  boost.Boost(p);
  k->part->SetMomentum(p);
  RoBoFin(k->left,rot,boost);
  RoBoFin(k->right,rot,boost);
  k->CheckMomentumConservation();
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
	msg_Error()<<METHOD<<"(..): Relative mass deviation "
		   <<(tnew/si->t-1.0)<<", "<<tnew<<" vs. "<<si->t
		   <<". Reset t."<<std::endl;
      si->t=tnew;
    }
    mo->part->SetMomentum(v_mo);
    si->part->SetMomentum(v_si);
    si->E2 = sqr(v_si[0]);
  }
  msg_Debugging()<<"}\n";
}

double Spacelike_Kinematics::BoostInCMS(Tree **const trees,Tree *const tree,
					const int tree1,
					Knot *const k1,Knot *const k2) 
{
  double dir((double)k1->dir);
  if (dir==0.0) THROW(fatal_error,"No knot history found");
  Vec4D cms(k1->part->Momentum()+k2->part->Momentum());
  m_boost=Poincare(cms);
  trees[tree1]->BoRo(m_boost);
  trees[1-tree1]->BoRo(m_boost);
  if (tree!=NULL) tree->BoRo(m_boost);
  if (dir>0.0) m_rot=Poincare(k1->part->Momentum(),Vec4D::ZVEC);
  else m_rot=Poincare(k2->part->Momentum(),Vec4D::ZVEC);
  
  trees[tree1]->BoRo(m_rot);
  trees[1-tree1]->BoRo(m_rot);
  if (tree!=NULL) tree->BoRo(m_rot);
  return dir;
}

double Spacelike_Kinematics::BoostFromCMS(Tree **const trees,Tree *const tree) 
{
  m_rot.Invert();
  trees[0]->BoRo(m_rot);
  trees[1]->BoRo(m_rot);
  if (tree!=NULL) tree->BoRo(m_rot);
  m_boost.Invert();
  trees[0]->BoRo(m_boost);
  trees[1]->BoRo(m_boost);
  if (tree!=NULL) tree->BoRo(m_boost);
  m_rot.Invert();
  m_boost.Invert();
  return 0.;
}

Vec4D Spacelike_Kinematics::BoostInLab(Tree **const trees,
				       Tree *const tree) 
{
  double x1(trees[0]->GetInitiator()->x), x2(trees[1]->GetInitiator()->x);
  // Only for massless initiators.
  Vec4D  lab(Vec4D(x1+x2,0.,0.,x2-x1));
  m_boost = Poincare(lab);
  trees[0]->BoRo(m_boost);
  trees[1]->BoRo(m_boost);
  if (tree!=NULL) tree->BoRo(m_boost);
  return trees[0]->GetRoot()->part->Momentum()+ 
    trees[1]->GetRoot()->part->Momentum();
}

bool Spacelike_Kinematics::ResetEnergies(Knot *const in) 
{
  double tt((in->stat!=3||in->shower>1)?in->t:in->tout);
  msg_Debugging()<<METHOD<<"("<<in->kn_no
		 <<"): test E = "<<sqrt(in->E2)<<" > m = "<<sqrt(tt)<<"\n";
  msg_Indent();
  if (in->E2 < tt) return 0;
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
					   Knot *const partner,bool diced) 
{
  if (active==NULL || partner==NULL || active->prev==NULL) {
    msg_Error()<<METHOD<<"("<<(active?(active->prev?active->prev->kn_no:-1):-1)
	       <<"->"<<(active?active->kn_no:-1)<<","
	       <<(partner?partner->kn_no:-1)<<"): Missing knot."<<std::endl;
    return 0.;
  }
  double t1(active->t), t2(partner->t), t3(diced?active->prev->t:active->prev->tout);
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
    msg_Debugging()<<"knots "<<active->prev->kn_no<<"->"<<active->kn_no<<","
		   <<active->prev->left->kn_no<<", maxt(3) = "
		   <<(maxt0+maxt1)<<" vs. "<<active->prev->left->tout
		   <<" <- t1 = "<<t1<<", t2 = "<<t2<<", t3 = "<<t3
		   <<", z = "<<active->z<<"\n";
    return maxt0+maxt1;
  }

  msg_Debugging()<<"knots "<<active->prev->kn_no<<"->"<<active->kn_no<<","
		 <<active->prev->left->kn_no<<", maxt(2) = "
		 <<(t1 + t3 + (np1*np3 - s1*s3)/(2.*t2))
		 <<" vs. "<<active->prev->left->tout
		 <<" <- t1 = "<<t1<<", t2 = "<<t2<<", t3 = "<<t3
		   <<", z = "<<active->z<<"\n";
  if (dabs(t2)>rpa.gen.Accu()) return  (t1 + t3 + (np1*np3 - s1*s3)/(2.*t2));
  msg_Debugging()<<"knots "<<active->prev->kn_no<<"->"<<active->kn_no<<","
		 <<active->prev->left->kn_no<<", maxt(1) = "
		 <<maxt0<<" vs. "<<active->prev->left->tout
		 <<" <- t1 = "<<t1<<", t2 = "<<t2<<", t3 = "<<t3
		 <<", z = "<<active->z<<"\n";
  return maxt0;
}

void Spacelike_Kinematics::ResetMomenta(Knot *const k, Tree *const tree, 
					const int &mode) 
{
  msg_Debugging()<<METHOD<<"(): "<<k->kn_no<<"\n";
  msg_Indent();
  if (mode==0) {
    Knot *mo(k->prev), *si(mo->left);
    // kill not-ME daughters of sister
    ResetMomenta(si,tree,1);
    // kill not-ME grand mother
    ResetMomenta(mo,tree,2);
    k->stat=1;
  }
  else if (mode==1) {
    if (k->part->Info()!='H' && k->decay==NULL) {
      k->left=NULL;
      k->right=NULL;
    }
    else {
      k->Restore(1);
      Vec4D mom(k->part->Momentum());
      if (k->decay!=NULL) {
	k->left=k->decay->left;
	k->right=k->decay->right;
	k->decay->left->prev=k;
	k->decay->right->prev=k;
      }
      tree->Restore(k);
      k->part->SetMomentum(mom);
      if (k->left) {
	if (k->left->part->Info()!='H' && k->left->decay==NULL) {
	  k->left=NULL;
	  k->right=NULL;
	}
	else {
	  ResetMomenta(k->left,tree,1);
	  ResetMomenta(k->right,tree,1);
	}
      }
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
