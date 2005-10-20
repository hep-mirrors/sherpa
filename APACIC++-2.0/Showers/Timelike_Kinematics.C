#include "Timelike_Kinematics.H"

#include "Knot.H"
#include "Tree.H"
#include "Poincare.H"
#include "Exception.H"

using namespace APACIC;
using namespace ATOOLS;

Timelike_Kinematics::Timelike_Kinematics(ATOOLS::Jet_Finder *const jf): 
  p_jf(jf),
  m_zscheme(1), m_anglescheme(0) {}

bool Timelike_Kinematics::CheckZRange(Knot * const  mo,
				      const Flavour *const d1flav,
				      const Flavour *const d2flav) const
{
  Knot * d1(mo->left), * d2(mo->right);
  if ((d1->stat==3) || (d2->stat==3)) return true; 

  double t(mo->t), t1(d1->t), t2(d2->t);

  if (t<t1+t2+2.*sqrt(t1*t2)) {
    if (d1->stat==0 && d2->stat==0)   return false;
    if (d1->stat==0 && d2->stat!=0) {
      d2->stat = 3;
      return true;
    }
    if (d1->stat!=0 && d2->stat==0) {
      d1->stat = 3;
      return true;
    }
    if (d1->t>d2->t) d1->stat = 3;
                else d2->stat = 3;
    return true;
  }

  double t01(d1->tout), t02(d2->tout);

  double z(GetZShift(mo->z,t,t1,t2,t01,t02)), 
    e12(z*z*mo->E2), e22((1.-z)*(1.-z)*mo->E2);
  bool do1(false), do2(false);
  if (d1->stat>0 && !CheckZRange(d1->z,e12,t1,sqr(d1flav[0].PSMass()),
				 sqr(d1flav[1].PSMass()))) do1=true;
  if (d2->stat>0 && !CheckZRange(d2->z,e22,t2,sqr(d2flav[0].PSMass()),
				 sqr(d2flav[1].PSMass()))) do2=true;

  if (!do1 && !do2) return true;
  else if (!do1 && do2) d2->stat=3;
  else if (do1 && !do2) d1->stat=3;
  else if (d1->t>d2->t) d1->stat=3;
                   else d2->stat=3;
  return false;
}

bool Timelike_Kinematics::
CheckZRange(const double &z, const double &E2, const double &t, 
	    const double &t1, const double &t2) const
{
  double x_p_mom(sqrt(1.-t/E2));   
  double mean_z(0.),delta_z(1.);              
  
  switch (m_zscheme) {
  case 0 : 
    mean_z  = 0.5;
    delta_z = 0.5*x_p_mom;
    break;
  case 1 :
    mean_z  = 0.5 *( 1. + (t1-t2)/t); 
    delta_z = 0.5 * x_p_mom * sqrt( sqr(t-t1-t2) - 4.*t1*t2)/t;
    break;
  case 2 :
    mean_z  = 0.5 * (1. + (t1-t2)/t); 
    delta_z = 0.5 * sqrt( sqr(t-t1-t2) - 4.*t1*t2)/t;    
    break;
  }
  
  if ((z<mean_z-delta_z) || (z>mean_z+delta_z)) {
    return false;
  }
  
  return true;
}

double Timelike_Kinematics::
GetZShift(const double &z, const double &t, 
	  const double &t1, const double &t2, 
	  const double &t01, const double &t02) const
{
  double lambda1(sqr(t-t1-t2)-4.*t1*t2); 
  if (lambda1<0) {
    msg_Tracking()<<"Timelike_Kinematics::CalcZShift(..): "
		  <<"kinematics does not fit : "
		  <<t<<" -> "<<t1<<"+"<<t2<<std::endl;
    return -1.;
  }
  switch (m_zscheme) {
  case 1:
    return ((2.*z-1.)*sqrt(lambda1) + (t+t1-t2))/(2.*t);
  case 0:
    double lambda0(sqr(t-t01-t02)-4.*t01*t02); 
    return  (z-(t+t01-t02)/(2.*t))*sqrt(lambda1/lambda0) + (t+t1-t2)/(2.*t);
  }
  return -1.0;
}

double Timelike_Kinematics::
GetZShift(const double &z, const double &t, 
	  const double &t1, const double &t2) const
{
  switch (m_zscheme) {
  case 0:
    return z;
  case 1:
    double lambda1 = sqr(t-t1-t2)-4.*t1*t2; 
    if (lambda1<0) {
      msg_Tracking()<<"Timelike_Kinematics::CalcZShift(..) "
		    <<"kinematics does not fit : "
		    <<t<<" -> "<<t1<<"+"<<t2<<std::endl;
      return -1.;
    }
    return ((2.*z-1.)*sqrt(lambda1) + (t+t1-t2))/(2.*t);
  }
  return -1.0;
}

int Timelike_Kinematics::Shuffle(Knot * const mo, const int first) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<first<<"): {\n";
  msg_Indent();
  int stat(1);
  if (first) stat=ShuffleMomenta(mo);
  else stat=ShuffleZ(mo);
  msg_Debugging()<<"z = "<<mo->z<<"\n";
  if (stat==1) stat=UpdateDaughters(mo);
  msg_Debugging()<<"}\n";
  return stat;
} 

int Timelike_Kinematics::ShuffleZ(Knot * const mo) const
{
  double t(mo->t);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout); 
  msg_Debugging()<<"t = "<<t<<", t_1 = "<<t1<<", t_2 = "<<t2<<"\n";
  if (t1+t2+2.*sqrt(t1*t2)-t>rpa.gen.Accu()) {
    msg_Debugging()<<METHOD<<"(..): Missing mass. m_a = "<<sqrt(t)
		   <<", m_b "<<sqrt(t1)<<", m_c = "<<sqrt(t2)<<std::endl;
    return 0; 
  }
  double t01(mo->left->tout), t02(mo->right->tout), z(mo->z);
  mo->z=GetZShift(z,t,t1,t2,t01,t02);
  if (!CheckKinematics(mo,0)) {
    mo->z = z;
    msg_Debugging()<<"failed: z = "<<mo->z<<", E_1 = "<<sqrt(mo->left->E2)
		   <<", E_2 = "<<sqrt(mo->right->E2)<<"\n";
    return 0;
  }
  mo->left->E2  = mo->z*mo->z*mo->E2;
  mo->right->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
  if (mo->z<=0.0 || mo->z>=1.0 ||
      mo->left->E2<mo->left->tout || mo->right->E2<mo->right->tout) return 0;
  msg_Debugging()<<"z = "<<mo->z<<", E_1 = "<<sqrt(mo->left->E2)
		 <<", E_2 = "<<sqrt(mo->right->E2)<<"\n";
  return 1;
}

int Timelike_Kinematics::ShuffleMomenta(Knot *const mo) const
{ 
  Knot *d1(mo->left), *d2(mo->right);
  double t1(d1->t), t2(d2->t);
  double ta(mo->part->Momentum().Abs2());
  if (dabs((mo->t-ta)/mo->t)>1.e-7) {
    msg.Error()<<METHOD<<"(..): Inconsistent masses. t = "<<mo->t
	       <<", p^2 = "<<ta<<std::endl;
  }
  double t(mo->t=ta);
  if (t1+t2+2.0*sqrt(t1*t2)-t>rpa.gen.Accu()) {
    msg_Debugging()<<"missing mass\n";
    return 0;
  }
  double r1(0.0), r2(0.0), z(mo->z);
  Vec4D p1(d1->part->Momentum()), p2(d2->part->Momentum());
  if (p1.Abs2()/d1->E2>rpa.gen.Accu() || p2.Abs2()/d2->E2>rpa.gen.Accu()) {
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
  if (dabs(mo->z-z) < rpa.gen.Accu()) {
    msg_Debugging()<<"shift unnecessary\n";
    // boost daughters
    BoostDaughters(mo);
    // update daughters 
    Tree::UpdateDaughters(mo);
    return 1;
  }
  if (!CheckKinematics(mo,1)) {
    msg_Debugging()<<"kinematics check failed "<<mo->z<<" "<<z<<std::endl;
    mo->z = z;
    return 0;
  }
  msg_Debugging()<<"z = "<<mo->z<<", p_1 = "<<d1->part->Momentum()
		 <<", p_2 = "<<d2->part->Momentum()<<"\n";
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
		 <<", p_2 = "<<d2->part->Momentum()<<"\n";
  // boost daughters
  BoostDaughters(mo);
  // update daughters 
  Tree::UpdateDaughters(mo);
  return 1;
}

int Timelike_Kinematics::UpdateDaughters(Knot *const mo,
					  const bool force) const
{
//   msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<force<<"): "
// 		 <<mo->part->Momentum()<<"\n";
//   msg_Indent();
  if (mo->left) {
    mo->left->E2=sqr(mo->z)*mo->E2;
    mo->right->E2=sqr((1.-mo->z))*mo->E2;
    int stat(1);
    if (mo->left->left!=NULL && mo->left->stat>0) stat=Shuffle(mo->left,0);
    if (stat==1 && 
	mo->right->left!=NULL && mo->right->stat>0) stat=Shuffle(mo->right,0);
    if (stat!=1) {
      msg_Debugging()<<METHOD<<"(..): shuffle failed\n";
      return stat;
    }
    DoSingleKinematics(mo,force);
    UpdateDaughters(mo->left);
    UpdateDaughters(mo->right);
  }
  return 1;
}

bool Timelike_Kinematics::CheckKinematics(Knot *const mo,
					  const int &first) const
{
  Knot *d1(mo->left), *d2(mo->right);
  if (d1==NULL || d2==NULL) return true;
  double E12(mo->z*mo->z*mo->E2), E22((1.0-mo->z)*(1.0-mo->z)*mo->E2);
  double t1(mo->left->stat!=3?mo->left->t:mo->left->tout); 
  double t2(mo->right->stat!=3?mo->right->t:mo->right->tout); 
  // timelike daughters 
  if (t1>E12 || t2>E22) {
    msg_Debugging()<<"timelike daughters\n";
    return false;
  }
  double p1p2(sqrt((E12-t1)*(E22-t2)));
  // triangular three momementum relation     
  if (mo->E2-mo->t-(E12-t1+E22-t2+2.0*p1p2)>rpa.gen.Accu()) {
    msg_Debugging()<<"three momentum\n";
    return false;
  }
  double costh((2.0*mo->z*(1.0-mo->z)*mo->E2-mo->t+t1+t2)/(2.*p1p2)); 
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
  if (!IsZero(mo->E2-mo->t)) {
    costh1=(2.0*mo->z*mo->E2-mo->t-t1+t2)/(2.0*sqrt((mo->E2-mo->t)*(E12-t1)));
  }
  if (dabs(costh1)>1.+rpa.gen.Accu()) {
    msg_Debugging()<<" |cos_1| = "<<costh1<<"\n";
    return false;
  }
  return true;
}

void Timelike_Kinematics::DoSingleKinematics(Knot * const mo,
					     const bool force) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<force<<"): "
		 <<mo->left->didkin<<"\n";
  if (!mo->left->didkin || force) {
    Vec4D p1(Vec4D(0.,0.,0.,0.)), p2(Vec4D(0.,0.,0.,0.));
    ConstructVectors(mo,p1,p2);
    mo->left->part->SetMomentum(p1);
    mo->right->part->SetMomentum(p2);
  }
}

bool Timelike_Kinematics::DoKinematics(Knot * const mo) const
{
  msg_Debugging()<<"Timelike_Kinematics::DoKinematics("
		 <<mo->kn_no<<"): {"<<std::endl;
  msg_Indent();
  if (!mo) return true;
  if (!mo->left) {
    if (mo->part->Info()==' ') {
      mo->part->SetStatus(1);
      mo->part->SetInfo('F');
    }
    msg_Debugging()<<"}\n";
    return true;
  }
  DoSingleKinematics(mo);
  mo->part->SetStatus(2);
  mo->left->didkin=true;
  mo->right->didkin=true;
  int error(0);
  if (!CheckVector(mo->part->Momentum()) || 
      !CheckVector(mo->left->part->Momentum()) || 
      !CheckVector(mo->right->part->Momentum())) error=1;
  static double accu(sqrt(rpa.gen.Accu()));
  Vec4D::SetAccu(accu);
  if (!(mo->left->part->Momentum()+mo->right->part->Momentum()==
	mo->part->Momentum())) error=2;
  Vec4D::ResetAccu();
  if (!(mo->left->part->Momentum()[2]<0) && 
      !(mo->left->part->Momentum()[2]>0)) error=3;
  if (error>0) {
    int op(msg.Error().precision(6));
    msg.Error()<<"Timelike_Kinematics::DoKinematics("<<mo->kn_no<<"): "
	       <<"Error "<<error<<": Momentum conservation violated.\n"
	       <<"   p      = "<<mo->part->Momentum()<<" -> "
	       <<mo->part->Momentum().Abs2()<<" vs. "<<mo->t
	       <<" ("<<mo->part->Flav()<<","<<mo->part->Info()<<")\n"
	       <<"   p_1    = "<<mo->left->part->Momentum()<<" -> "
	       <<mo->left->part->Momentum().Abs2()<<" vs. "
	       <<mo->left->t<<" ("<<mo->left->part->Flav()<<","
	       <<mo->left->part->Info()<<")\n"
	       <<"   p_2    = "<<mo->right->part->Momentum()<<" -> "
	       <<mo->right->part->Momentum().Abs2()<<" vs. "
	       <<mo->right->t<<" ("<<mo->right->part->Flav()<<","
	       <<mo->right->part->Info()<<")\n"
	       <<"   p_miss = "<<(mo->part->Momentum()-
				  mo->left->part->Momentum()-
				  mo->right->part->Momentum())<<std::endl;
    msg.Error().precision(op);
    return false;
  }
  if (!DoKinematics(mo->left)) return false;
  if (!DoKinematics(mo->right)) return false;
  BoostDaughters(mo);
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
  double p1(sqrt(E12-t1)), p2(sqrt(E22-t2));
//   PRINT_INFO(t1<<" "<<t2<<" "<<p1<<" "<<p2<<" "<<p);
  Vec3D n1,n2;
  ConstructDreiBein(mo,n1,n2);
  double phi(mo->phi + mo->polinfo.Angle()), bph(cos(phi)), cph(-sin(phi));
  Vec3D es(cph*n1 + bph*n2);
  double cth1((p*p-p2*p2+p1*p1)/(2.*p*p1)), sth1(sqrt(1.-sqr(cth1)));
  double cth2((p*p+p2*p2-p1*p1)/(2.*p*p2)), sth2(sqrt(1.-sqr(cth2)));
  mo->costh   = cth1*cth2-sth1*sth2;
  Vec3D nm(mo->part->Momentum());
  nm=1.0/nm.Abs()*nm;
  p1vec = Vec4D(sqrt(E12),p1*(cth1*nm - sth1*es));
  p2vec = Vec4D(sqrt(E22),p2*(cth2*nm + sth2*es));
  if (p1vec.Nan() || p2vec.Nan()) {
    msg.Error()<<METHOD<<"("<<mo->kn_no<<"): Error."<<std::endl;
    msg_Debugging()<<"mo = "<<*mo;
    msg_Debugging()<<"d1 = "<<*mo->left;
    msg_Debugging()<<"d2 = "<<*mo->right;
    msg_Debugging()<<nm<<" "<<es<<" "<<phi<<" "<<cth1<<" "<<cth2<<std::endl;
  }
}

void Timelike_Kinematics::
ConstructDreiBein(Knot *const mo,Vec3D &n1,Vec3D &n2) const
{
  if (mo->prev==NULL) {
    n1 = Vec3D(0.,0.,1.);
    n2 = Vec3D(0.,1.,0.);
    return;
  }
  Knot * au(mo->prev->left);
  int sign(0), mode(1);
  if (mo==au) {
    au      = mo->prev->right;
    sign    = 1;
    mode    = 3;
  }
  Vec3D na(au->part->Momentum()); // aunt
  Vec3D nm(mo->part->Momentum()); // mother
  n1 = cross(na,nm);
  double n1abs(n1.Abs());
  if (n1abs<1.e-5) {
    n1   = cross(Vec3D(0.,0.,1.),nm);
    mode = mode|4;
  }
  if (n1abs<1.e-5) {
    n1   = cross(Vec3D(0.,1.,0.),nm);
    mode = mode|8;
  }  
  if (sign) n1 = -1.*n1; 
  n2 = cross(nm,n1);
  n1 = n1/n1.Abs();
  n2 = n2/n2.Abs();
}

bool Timelike_Kinematics::CheckVector(const Vec4D &mom) const 
{
  double abs2(mom.Abs2());
  if (abs2>0 && abs2<0) return false;
  if (mom[0]<0) return false;
  return true;
}
 
void Timelike_Kinematics::BoostDaughters(Knot * const mo) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):\n";
  msg_Indent();
  Knot *d1(mo->left), *d2(mo->right);
  static double accu(sqrt(rpa.gen.Accu()));
  Vec4D::SetAccu(accu);
  if (d1->left && d1->left->part->Momentum()!=Vec4D()) {
    Poincare cms(d1->left->part->Momentum()+d1->right->part->Momentum());
    Poincare cmsp(d1->part->Momentum());
    cmsp.Invert();
    Vec4D tp(d1->left->part->Momentum()+d1->right->part->Momentum());
    cms.Boost(tp);
    cmsp.Boost(tp);
    Tree::BoRo(cms,d1->left);
    Tree::BoRo(cms,d1->right);
    Tree::BoRo(cmsp,d1->left);
    Tree::BoRo(cmsp,d1->right);
    if (!(d1->part->Momentum()==
	  d1->left->part->Momentum()+d1->right->part->Momentum())) {
      msg.Error()<<METHOD<<"(..): Four momentum not conserved."
		 <<"\n p_miss = "<<(d1->part->Momentum()-
				    d1->left->part->Momentum()-
				    d1->right->part->Momentum())
		 <<"\n p      = "<<d1->part->Momentum()
		 <<"\n p_1    = "<<d1->left->part->Momentum()
		 <<"\n p_2    = "<<d1->right->part->Momentum()
		 <<"."<<std::endl;
    }
  }
  if (d2->left && d2->left->part->Momentum()!=Vec4D()) {
    Poincare cms(d2->left->part->Momentum()+d2->right->part->Momentum());
    Poincare cmsp(d2->part->Momentum());
    cmsp.Invert();
    Tree::BoRo(cms,d2->left);
    Tree::BoRo(cms,d2->right);
    Tree::BoRo(cmsp,d2->left);
    Tree::BoRo(cmsp,d2->right);
    if (!(d2->part->Momentum()==
	  d2->left->part->Momentum()+d2->right->part->Momentum())) {
      msg.Error()<<METHOD<<"(..): Four momentum not conserved."
		 <<"\n p_miss = "<<(d2->part->Momentum()-
				    d2->left->part->Momentum()-
				    d2->right->part->Momentum())
		 <<"\n p      = "<<d2->part->Momentum()
		 <<"\n p_1    = "<<d2->left->part->Momentum()
		 <<"\n p_2    = "<<d2->right->part->Momentum()
		 <<"."<<std::endl;
    }
  }
  Vec4D::ResetAccu();
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
  double costh((2.0*z*(1.0-z)*E2-ta+tb+tc)/(2.0*sqrt((E12-tb)*(E22-tc))));
  if (dabs(costh)>1.0) {
//     msg.Error()<<METHOD<<"(..): cos_{theta} = "<<costh<<std::endl;
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
//     msg.Error()<<METHOD<<"(..): cos_{theta} = "<<costh<<std::endl;
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
//   if (zlc<0.0 || zlc>1.0) {
//     msg.Error()<<METHOD<<"(..): z = "<<zlc<<std::endl;
//   }
  return zlc;
}

bool Timelike_Kinematics::
ArrangeColourPartners(Particle *const aup,Knot *const d1,Knot *const d2) const
{
  if (!aup || !d1 || !d2)  return false;
  int jfmode(p_jf->Type());
  p_jf->SetType(1);
  bool left(p_jf->MTij2(aup->Momentum(),d1->part->Momentum())<
	    p_jf->MTij2(aup->Momentum(),d2->part->Momentum()));
  p_jf->SetType(jfmode);
  if (left) return false;
  return true;
}

