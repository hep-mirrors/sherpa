#include "Timelike_Kinematics.H"

#include "Knot.H"
#include "Tree.H"
#include "Poincare.H"
#include "Exception.H"

using namespace APACIC;
using namespace ATOOLS;

Timelike_Kinematics::Timelike_Kinematics(ATOOLS::Jet_Finder *const jf): 
  p_jf(jf),
  m_zscheme(1) {}

double Timelike_Kinematics::GetZ(const double &zp, const double &t, 
				 const double &t1, const double &t2, 
				 const double &t01, const double &t02) const
{
  double lambda1(sqr(t-t1-t2)-4.0*t1*t2); 
  if (lambda1<0.0) {
    msg_Tracking()<<METHOD<<"(..): Bad kinematics. t = "
		  <<t<<", t_1 = "<<t1<<", t_2 = "<<t2<<std::endl;
    return -1.0;
  }
  switch (m_zscheme) {
  case 1:
    return ((2.0*zp-1.0)*sqrt(lambda1)+(t+t1-t2))/(2.0*t);
  case 0:
    double lambda0(sqr(t-t01-t02)-4.0*t01*t02); 
    return (zp-(t+t01-t02)/(2.0*t))*sqrt(lambda1/lambda0)+(t+t1-t2)/(2.0*t);
  }
  return -1.0;
}

int Timelike_Kinematics::Shuffle(Knot * const mo, const int first) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<first<<"): {\n";
  msg_Indent();
  int stat(first&&(mo->shower<2||mo->stat==0)?ShuffleMomenta(mo):ShuffleZ(mo));
  msg_Debugging()<<"z = "<<mo->z<<", stat = "<<stat<<"\n";
  if (stat==1) stat=UpdateDaughters(mo);
  msg_Debugging()<<"}\n";
  return stat;
} 

int Timelike_Kinematics::UpdateDaughters(Knot *const mo,
					 const bool force) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<force<<")\n";
  if (mo->left) {
    int stat(1);
    if (mo->left->left!=NULL && mo->left->stat!=3)
      stat=Shuffle(mo->left,mo->left->left->part->Info()=='H' ||
		   mo->left->right->part->Info()=='H');
    if (stat==1 && 
	mo->right->left!=NULL && mo->right->stat!=3) 
      stat=Shuffle(mo->right,mo->right->left->part->Info()=='H' ||
		   mo->right->right->part->Info()=='H');
    if (stat!=1) {
      msg_Debugging()<<METHOD<<"(..): shuffle failed\n";
      mo->z=mo->zs;
      mo->left->E2=sqr(mo->z)*mo->E2;
      mo->right->E2=sqr((1.-mo->z))*mo->E2;
      return stat;
    }
  }
  return 1;
}

int Timelike_Kinematics::GeneratePSMasses(Knot *const mo) const
{
  if (mo->oc[0]<0) {
    mo->oc[0]=mo->part->GetFlow(1);
    mo->oc[1]=mo->part->GetFlow(2);
  }
  if (mo->left==NULL) return 1;
  msg_Debugging()<<METHOD<<"(): knot "<<mo->kn_no<<"\n";
  mo->E2=sqr(mo->part->Momentum()[0]);
  msg_Indent();
  int res(1);
  if ((res=GeneratePSMasses(mo->left))!=1) return res;
  if ((res=GeneratePSMasses(mo->right))!=1) return res;
  mo->z=mo->left->part->Momentum()[0]/mo->part->Momentum()[0];
  if (mo->left->left==NULL || mo->right->left==NULL) 
    res=ShuffleMomenta(mo,true,true);
  mo->CheckMomentumConservation(); 
  if (!res) {
    msg_Tracking()<<METHOD<<"(): Invalid kinematics in splitting "
		  <<mo<<"->("<<mo->left<<","<<mo->right<<") {\n\n"
		  <<*mo<<*mo->left<<*mo->right<<"}"<<std::endl;
  }
  return res;
}

int Timelike_Kinematics::ShuffleZ(Knot * const mo,const bool update) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<")\n";
  msg_Indent();
  double t(mo->t);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout); 
  double t01(mo->left->tout), t02(mo->right->tout);
  if (mo->shower==3) t=mo->tmo;
  if (mo->left->shower==3) t01=t1=mo->left->tmo;
  if (mo->right->shower==3) t02=t2=mo->right->tmo;
  msg_Debugging()<<"t = "<<t<<", t_1 = "<<t1<<", t_2 = "<<t2<<"\n";
  if (t1+t2+2.0*sqrt(t1*t2)-t>rpa.gen.Accu()) {
    msg_Debugging()<<METHOD<<"(..): Missing mass. m_a = "<<sqrt(t)
		   <<", m_b "<<sqrt(t1)<<", m_c = "<<sqrt(t2)<<std::endl;
    return 0; 
  }
  msg_Debugging()<<"z = "<<mo->zs<<", t0_1 = "<<t01<<", t0_2 = "<<t02<<"\n";
  mo->z=GetZ(mo->zs,t,t1,t2,t01,t02);
  if (!CheckKinematics(mo,0)) {
    mo->z=mo->zs;
    msg_Debugging()<<"failed: z = "<<mo->z<<", E_1 = "<<sqrt(mo->left->E2)
		   <<", E_2 = "<<sqrt(mo->right->E2)<<"\n";
    return 0;
  }
  mo->left->E2=sqr(mo->z)*mo->E2;
  mo->right->E2=sqr(1.0-mo->z)*mo->E2;
  static double accu(sqrt(rpa.gen.Accu())); 
  if (mo->z<accu || 1.0-mo->z<accu || 
      mo->left->E2<mo->left->tout || mo->right->E2<mo->right->tout) return 0;
  msg_Debugging()<<"z = "<<mo->z<<", 1-z = "<<1.0-mo->z
		 <<", E_1 = "<<sqrt(mo->left->E2)<<", E_2 = "
		 <<sqrt(mo->right->E2)<<", E = "<<sqrt(mo->E2)<<"\n";
  if (mo->left->part->Info()!='H') mo->right->didkin=mo->left->didkin=false;
  DoSingleKinematics(mo,true);
  if (update) {
    if (!BoostDaughters(mo)) return 0;
    Tree::UpdateDaughters(mo);
  }
  mo->CheckMomentumConservation();
  return 1;
}

int Timelike_Kinematics::ShuffleMomenta
(Knot *const mo,const bool update,const bool gpm) const
{ 
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"->"
		 <<mo->left->kn_no<<","<<mo->right->kn_no<<")\n";
  msg_Indent();
  Knot *d1(mo->left), *d2(mo->right);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout); 
  if (mo->left->shower==3) t1=mo->left->tmo;
  if (mo->right->shower==3) t2=mo->right->tmo;
  double t(mo->part->Momentum().Abs2());
  msg_Debugging()<<"t = "<<t<<", t_1 = "<<t1<<", t_2 = "<<t2<<"\n";
  if (dabs((mo->t-t)/mo->t)>1.e-7 && mo->shower!=2) {
    //msg_Error()<<METHOD<<"(..): Inconsistent masses. t = "<<mo->t<<", p^2 = "<<t<<std::endl;
    mo->t=t;
  }
  if (t1+t2+2.0*sqrt(t1*t2)-t>rpa.gen.Accu()) {
    msg_Debugging()<<"missing mass "<<sqrt(t1+t2+2.0*sqrt(t1*t2))
		   <<" vs. "<<sqrt(t)<<"\n";
    return 0;
  }
  double r1(0.0), r2(0.0), z(mo->z);
  Vec4D p1(d1->part->Momentum()), p2(d2->part->Momentum());
  if ((t1+p1.Abs2())/d1->E2>rpa.gen.Accu() || 
      (t2+p2.Abs2())/d2->E2>rpa.gen.Accu()) {
    double t1n(p1.Abs2()), t2n(p2.Abs2());
    double A(((t2-t2n)-(t1-t1n))/(t+t1n-t2n));
    double B((t+t2n-t1n)/(t+t1n-t2n));
    double C(t-t1n-t2n);
    double D((2.0*t2n-2.0*A*B*t1n+(A-B)*C)/(2.0*(t2n+B*B*t1n-B*C)));
    double E((t2n-t2+A*A*t1n+A*C)/(t2n+B*B*t1n-B*C));
    r2=D-sqrt(D*D-E);
    r1=A+r2*B;
    Vec4D p1a((1.0-r1)*p1+r2*p2);
    Vec4D p(mo->part->Momentum());
    mo->z=p1a[0]/p[0];
  }
  else {
    double lambda(sqrt(sqr(t-t1-t2)-4.0*t1*t2)); 
    r1=(t+t2-t1-lambda)/(2.0*t);
    r2=(t-t2+t1-lambda)/(2.0*t);
    mo->z=z-r1*z+r2*(1.0-z);
  } 
  if (!CheckKinematics(mo,1,gpm)) {
    msg_Debugging()<<"kinematics check failed "<<mo->z<<" "<<z<<std::endl;
    mo->z = z;
    return 0;
  }
  msg_Debugging()<<"z = "<<mo->z<<", p_1 = "<<d1->part->Momentum()
		 <<" "<<d1->part->Momentum().Abs2()
		 <<", p_2 = "<<d2->part->Momentum()
		 <<" "<<d2->part->Momentum().Abs2()<<"\n";
  msg_Debugging()<<"p-p_1-p_2, old = "<<(mo->part->Momentum()-
					 d1->part->Momentum()-
					 d2->part->Momentum())<<"\n";
  d1->part->SetMomentum((1.0-r1)*p1+r2*p2);
  d2->part->SetMomentum((1.0-r2)*p2+r1*p1);
  d1->E2=mo->z*mo->z*mo->E2;
  d2->E2=(1.0-mo->z)*(1.0-mo->z)*mo->E2;
  msg_Debugging()<<"           new = "<<(mo->part->Momentum()-
					 d1->part->Momentum()-
					 d2->part->Momentum())<<"\n";
  msg_Debugging()<<"z = "<<mo->z<<", p_1 = "<<d1->part->Momentum()
		 <<" "<<d1->part->Momentum().Abs2()
		 <<", p_2 = "<<d2->part->Momentum()<<" "
		 <<d2->part->Momentum().Abs2()<<"\n";
  if (update) {
    if (!BoostDaughters(mo)) {
      mo->z=mo->zs;
      d1->part->SetMomentum(p1);
      d2->part->SetMomentum(p2);
      d1->E2=mo->z*mo->z*mo->E2;
      d2->E2=(1.0-mo->z)*(1.0-mo->z)*mo->E2;
      return 0;
    }
    Tree::UpdateDaughters(mo);
  }
  return 1;
}

bool Timelike_Kinematics::CheckKinematics
(Knot *const mo,const int &first,const bool gpm) const
{
  Knot *d1(mo->left), *d2(mo->right);
  if (d1==NULL || d2==NULL) return true;
  double E12(mo->z*mo->z*mo->E2), E22((1.0-mo->z)*(1.0-mo->z)*mo->E2);
  double t(mo->t);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout); 
  if (mo->shower==3) t=mo->tmo;
  if (d1->shower==3) t1=d1->tmo;
  if (d2->shower==3) t2=d2->tmo;
  if (gpm) {
    t=mo->part->Momentum().Abs2();
    if (d1->left!=NULL) t1=d1->part->Momentum().Abs2();
    if (d2->left!=NULL) t2=d2->part->Momentum().Abs2();
  }
  msg_Debugging()<<"t = "<<t<<" ("<<mo->t<<"), t_1 = "<<t1
		 <<" ("<<d1->t<<"), t_2 = "<<t2<<" ("<<d2->t<<")\n";
  msg_Debugging()<<"E = "<<sqrt(mo->E2)<<", E_1 = "<<sqrt(E12)
		 <<", E_2 = "<<sqrt(E22)<<"\n";
  // timelike daughters 
  if (t1>E12 || t2>E22) {
    msg_Debugging()<<"timelike daughters\n";
    return false;
  }
  double p1p2(sqrt((E12-t1)*(E22-t2)));
  // triangular three momementum relation     
  if (mo->E2-t-(E12-t1+E22-t2+2.0*p1p2)>rpa.gen.Accu()) {
    msg_Debugging()<<"three momentum\n";
    return false;
  }
  double costh((2.0*mo->z*(1.0-mo->z)*mo->E2-t+t1+t2)/(2.0*p1p2)); 
  // physical opening angle
  if (dabs(costh)>1.0 && !first) {
    msg_Debugging()<<"|cos| > 1\n";
    return false;
  }
  if (costh>1.0) costh=1.0; 
  if (costh<-1.0) costh=-1.0; 
  mo->costh=costh;
  // physical deflection angle
  double costh1(-1.0);
  if (!IsZero(mo->E2-t)) {
    costh1=(2.0*mo->z*mo->E2-t-t1+t2)/(2.0*sqrt((mo->E2-t)*(E12-t1)));
  }
  if (dabs(costh1)>1.0+rpa.gen.Accu()) {
    msg_Debugging()<<" |cos_1| = "<<costh1<<"\n";
    return false;
  }
  return true;
}

bool Timelike_Kinematics::DoSingleKinematics
(Knot * const mo,const bool force) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<force<<"): "
		 <<mo->left->didkin<<"\n";
  msg_Indent();
  if (!mo->left->didkin || !mo->right->didkin || force) {
    Vec4D p1, p2;
    ConstructVectors(mo,p1,p2);
    mo->left->part->SetMomentum(p1);
    mo->right->part->SetMomentum(p2);
    mo->left->didkin=true;
    mo->right->didkin=true;
  }
  if (!mo->CheckMomentumConservation(true)) return false;
  if (!CheckVector(mo->part->Momentum()) || 
      !CheckVector(mo->left->part->Momentum()) || 
      !CheckVector(mo->right->part->Momentum())) {
    msg_Error()<<METHOD<<"(): Constructed negative energy momentum."
	       <<std::endl;
    return false;
  }
  return true;
}

bool Timelike_Kinematics::DoKinematics(Knot * const mo) const
{
  msg_Debugging()<<"Timelike_Kinematics::DoKinematics("
		 <<mo->kn_no<<"): {"<<std::endl;
  msg_Indent();
  if (!mo) return true;
  if (!mo->left) {
    if (mo->part->Info()==' ') {
      mo->part->SetStatus(part_status::active);
      mo->part->SetInfo('F');
    }
    msg_Debugging()<<"}\n";
    return true;
  }
  if (!DoSingleKinematics(mo)) return false;
  mo->part->SetStatus(part_status::decayed);
  if (!DoKinematics(mo->left)) return false;
  if (!DoKinematics(mo->right)) return false;
  msg_Debugging()<<"}\n";
  return true;
}


void Timelike_Kinematics::
ConstructVectors(Knot *const mo,Vec4D &p1vec,Vec4D &p2vec) const
{
  double p(sqrt(mo->E2-mo->t));
  double E12(mo->left->E2), E22(mo->right->E2);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout);
  if (mo->shower==3) p=sqrt(mo->E2-mo->tmo);
  if (mo->left->shower==3) t1=mo->left->tmo;
  if (mo->right->shower==3) t2=mo->right->tmo;
  msg_Debugging()<<"construct "<<mo->kn_no<<", t_1 = "<<t1<<"("
		 <<mo->left->kn_no<<","<<mo->left->part->Info()
		 <<") <- "<<mo->left->t<<", t_2 = "<<t2<<"("<<mo->right->kn_no
		 <<","<<mo->right->part->Info()<<") <- "<<mo->right->t<<"\n";
  if (mo->left->part->Info()=='H' && mo->right->part->Info()=='H') {
    Vec4D p(mo->part->Momentum());
    Vec4D p1(mo->left->part->Momentum()), p2(mo->right->part->Momentum());
    double ap(p.PSpat());
    double kt(sqrt(GetRelativeKT2(mo->z,mo->E2,mo->t,t1,t2)));
    Vec4D l(0.0,1.0/ap*(Vec3D)p), t(0.0,(Vec3D)(p1-(p1*l)*l));
    t=1.0/t.PSpat()*t;
    if (p1.PSpat2()>p2.PSpat2()) {
      double lf(sqrt(E12-kt*kt));
      p1vec=kt*t+lf*l;
      p2vec=(-kt)*t+(ap-lf)*l;
    }
    else {
      double lf(sqrt(E22-kt*kt));
      p2vec=(-kt)*t+lf*l;
      p1vec=kt*t+(ap-lf)*l;
    }
    p1vec[0]=sqrt(E12);
    p1vec[1]=sqrt(E22);
    msg_Debugging()<<"old p_1 = "<<mo->left->part->Momentum()
		   <<", p_2 = "<<mo->right->part->Momentum()<<"\n";
    msg_Debugging()<<"new p_1 = "<<p1vec<<", p_2 = "<<p2vec<<"\n";
    if (p!=p1vec+p2vec) {
      msg_Error()<<METHOD<<"(..): Four momentum not conserved.\n"
		 <<"  p_miss  = "<<(p1vec+p2vec-p)<<"\n"
		 <<"  p_old   = "<<(p1+p2)<<" "<<(p1+p2).Abs2()<<"\n"
		 <<"  p_new   = "<<(p1vec+p2vec)
		 <<" "<<(p1vec+p2vec).Abs2()<<std::endl;
    }
  }
  else {
    Vec3D n1,n2;
    ConstructDreiBein(mo,n1,n2);
    double p1(sqrt(E12-t1)), p2(sqrt(E22-t2));
    msg_Debugging()<<"construct, p_1 = "<<p1<<", p_2 = "<<p2<<", p = "<<p<<"\n";
    double phi(mo->phi + mo->polinfo.Angle()), bph(cos(phi)), cph(-sin(phi));
    Vec3D es(cph*n1 + bph*n2);
    double cth1((p*p-p2*p2+p1*p1)/(2.*p*p1)), sth1(sqrt(1.-sqr(cth1)));
    if (!(sth1>0.0) && IsZero(cth1-1.0)) sth1=0.0;
    double cth2((p*p+p2*p2-p1*p1)/(2.*p*p2)), sth2(sqrt(1.-sqr(cth2)));
    if (!(sth2>0.0) && IsZero(cth2-1.0)) sth2=0.0;
    mo->costh=cth1*cth2-sth1*sth2;
    Vec3D nm(mo->part->Momentum());
    nm=1.0/nm.Abs()*nm;
    p1vec=Vec4D(sqrt(E12),p1*(cth1*nm - sth1*es));
    p2vec=Vec4D(sqrt(E22),p2*(cth2*nm + sth2*es));
    if (p1vec.Nan() || p2vec.Nan()) {
      msg_Error()<<METHOD<<"("<<mo->kn_no<<"): Error."<<std::endl
		 <<"mo = "<<*mo<<"d1 = "<<*mo->left<<"d2 = "<<*mo->right
		 <<"n = "<<nm<<", e = "<<es<<"\nphi = "
		 <<phi<<" <- "<<mo->phi<<" ("<<mo->polinfo.Angle()<<")"
		 <<"sth_1 = "<<sth1<<" <- "<<cth1-1.0<<", sth_2 = "
		 <<sth2<<" <- "<<cth2-1.0<<std::endl;
    }
  }
}

void Timelike_Kinematics::
ConstructDreiBein(Knot *const mo,Vec3D &n1,Vec3D &n2) const
{
  if (mo->prev==NULL) {
    n1=Vec3D(0.,0.,1.);
    n2=Vec3D(0.,1.,0.);
    return;
  }
  Knot *au(mo->prev->left);
  int sign(0), mode(1);
  if (mo==au) {
    au=mo->prev->right;
    sign=1;
    mode=3;
  }
  Vec3D na(au->part->Momentum()); // aunt
  Vec3D nm(mo->part->Momentum()); // mother
  n1=cross(na,nm);
  double n1abs(n1.Abs());
  if (n1abs<1.e-5) {
    n1=cross(Vec3D(0.,0.,1.),nm);
    mode=mode|4;
  }
  if (n1abs<1.e-5) {
    n1=cross(Vec3D(0.,1.,0.),nm);
    mode=mode|8;
  }  
  if (sign) n1=-1.*n1; 
  n2=cross(nm,n1);
  n1=n1/n1.Abs();
  n2=n2/n2.Abs();
}

bool Timelike_Kinematics::BoostDaughter(Knot * const d) const
{
  msg_Debugging()<<METHOD<<"("<<d->kn_no<<"): "
		 <<d->kn_no<<"->("<<(d->left?d->left->kn_no:0)
		 <<","<<(d->right?d->right->kn_no:0)<<")\n";
  msg_Indent();
  if (d->left && d->left->shower==3) {
    if (!ShuffleZ(d,false)) return false;
    if (!ReconstructDaughter(d->left)) return false;
  }
  if (d->left && d->left->part->Momentum()!=Vec4D()) {
    Poincare cms(d->left->part->Momentum()+d->right->part->Momentum());
    Poincare cmsp(d->part->Momentum());
    cmsp.Invert();
    Tree::BoRo(cms,d->left);
    Tree::BoRo(cms,d->right);
    Tree::BoRo(cmsp,d->left);
    Tree::BoRo(cmsp,d->right);
  }
  d->CheckMomentumConservation(true);
  return true;
}

bool Timelike_Kinematics::ReconstructDaughter(Knot * const d) const
{
  Knot *d1(d->left), *d2(d->right);
  msg_Debugging()<<METHOD<<"("<<d->kn_no<<"): "
		 <<d->kn_no<<"->("<<(d1?d1->kn_no:0)
		 <<","<<(d2?d2->kn_no:0)<<")\n";
  msg_Indent();
  Vec4D op1(d1->part->Momentum()), op2(d2->part->Momentum());
  Poincare cms(op1+op2);
  Poincare cmsp(d->part->Momentum());
  cmsp.Invert();
  msg_Debugging()<<"lab: p1 = "<<d1->part->Momentum()
		 <<", p2 = "<<d2->part->Momentum()<<"\n";
  Tree::BoRo(cms,d1);
  Tree::BoRo(cms,d2);
  msg_Debugging()<<"cms: p1 = "<<d1->part->Momentum()
		 <<", p2 = "<<d2->part->Momentum()<<"\n";
  Vec3D n(d1->part->Momentum());
  n=1.0/n.Abs()*n;
  double s1(d1->part->Momentum().Abs2()), s2(d2->part->Momentum().Abs2());
  double s(d->t), eps(0.5*(s+s1-s2)/s), pi2(eps*eps-s1/s);
  if (pi2<0.0) {
    msg_Error()<<METHOD<<"(): Cannot construct new cms."<<std::endl;
    abort();
  }
  double z(eps+sqrt(pi2)), Q(sqrt(d->t));
  msg_Debugging()<<"now: n = "<<n<<" "<<n.Abs()<<", z = "<<z
		 <<", Q = "<<Q<<", s1 = "<<s1<<", s2 = "<<s2<<"\n";
  double p1p(z*Q), p1m(s1/p1p);
  Vec4D p1(0.5*(p1p+p1m),0.5*(p1p-p1m)*n);
  d1->part->SetMomentum(p1);
  d2->part->SetMomentum(Vec4D(Q,0.0,0.0,0.0)-p1);
  msg_Debugging()<<"new cms: p1 = "<<d1->part->Momentum()
		 <<" "<<d1->part->Momentum().Abs2()
		 <<", p2 = "<<d2->part->Momentum()
		 <<" "<<d2->part->Momentum().Abs2()<<"\n";
  if (!BoostDaughter(d1) || !BoostDaughter(d2)) {
    d1->part->SetMomentum(op1);
    d2->part->SetMomentum(op2);
    return false;
  }
  d1->CheckMomentumConservation();
  d2->CheckMomentumConservation();
  Tree::BoRo(cmsp,d1);
  Tree::BoRo(cmsp,d2);
  msg_Debugging()<<"new lab: p1 = "<<d1->part->Momentum()
		 <<", p2 = "<<d2->part->Momentum()<<"\n";
  d->CheckMomentumConservation();
  return true;
}

bool Timelike_Kinematics::BoostDaughters(Knot * const mo) const
{
  Knot *d1(mo->left), *d2(mo->right);
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"): d1 = "
		 <<d1->kn_no<<", d2 = "<<d2->kn_no<<"\n";
  msg_Indent();
  if (!BoostDaughter(d1)) return false;
  if (!BoostDaughter(d2)) return false;
  mo->CheckMomentumConservation();
  return true;
}

double Timelike_Kinematics::GetOpeningAngle(Knot *const knot) const
{
  return GetOpeningAngle(knot->z,knot->E2,knot->t,
			 knot->left->t,knot->right->t);
}

double Timelike_Kinematics::
GetOpeningAngle(const double &z,const double &E2,const double &ta,
		const double &tb,const double &tc) const
{
  // opening angle of daughter partons in current frame
  double E12(z*z*E2), E22((1.0-z)*(1.0-z)*E2);
  if (E12<tb || E22<tc) return 0.0;
  //exact: 
  double costh((2.0*z*(1.0-z)*E2-ta+tb+tc)/(2.0*sqrt((E12-tb)*(E22-tc))));
  //approx: 
  //double costh((2.0*z*(1.0-z)*E2-ta)/(2.0*sqrt(E12*E22)));
  if (dabs(costh)>1.0) {
    return M_PI;
  }
  return acos(costh);
}

double Timelike_Kinematics::GetDeflectionAngle(Knot *const knot) const
{
  return GetDeflectionAngle(knot->z,knot->E2,knot->t,
			 knot->left->t,knot->right->t);
}

double Timelike_Kinematics::
GetDeflectionAngle(const double &z,const double &E2,const double &ta,
		const double &tb,const double &tc) const
{
  // deflection angle of first daughter parton in current frame
  double costh((2.0*z*E2-ta-tb+tc)/(2.0*sqrt((E2-ta)*(z*z*E2-tb))));
  if (dabs(costh)>1.0) {
    return M_PI;
  }
  return acos(costh);
}

double Timelike_Kinematics::
GetRelativeKT2(const double &z,const double &E2,
	       const double &ta,const double &tb,const double &tc) const
{
  // kt2 of daughter partons w.r.t. mother in light cone kinematics
  double zlc(LightConeZ(z,E2,ta,tb,tc));
  return zlc*(1.0-zlc)*ta-(1.0-zlc)*tb-zlc*tc;
}

double Timelike_Kinematics::LightConeZ(Knot *const knot) const
{
  return LightConeZ(knot->z,knot->E2,knot->t,knot->left->t,knot->right->t);
}

double Timelike_Kinematics::LightConeZ(const double &z,const double &E2,
				       const double &ta,
				       const double &tb,const double &tc) const
{
  // light cone momentum fraction of first daughter in light cone kinematics
  double pph(1.0+sqrt(1.0-ta/E2));
  double zlc(((ta+tb-tc)-2.0*z*E2*pph)/(ta-pph*pph*E2));
  return zlc;
}

double Timelike_Kinematics::EnergyZ(const double &zlc,const double &E2,
				    const double &ta,
				    const double &tb,const double &tc) const
{
  // energy fraction of first daughter
  double pph(1.0+sqrt(1.0-ta/E2));
  double z(((ta+tb-tc)-zlc*(ta-pph*pph*E2))/(2.0*E2*pph));
  return z;
}

bool Timelike_Kinematics::
ArrangeColourPartners(Particle *const aup,Knot *const d1,Knot *const d2) const
{
  if (!aup || !d1 || !d2)  return false;
  int jfmode(p_jf->Type());
  p_jf->SetType(1);
  bool left(p_jf->MTij2(aup->Momentum(),d1->part->Momentum(),
			aup->Flav().Mass(),d1->part->Flav().Mass())<
	    p_jf->MTij2(aup->Momentum(),d2->part->Momentum(),
			aup->Flav().Mass(),d2->part->Flav().Mass()));
  p_jf->SetType(jfmode);
  if (left) return false;
  return true;
}

