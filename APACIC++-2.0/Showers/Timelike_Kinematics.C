#include "Timelike_Kinematics.H"

#include "Run_Parameter.H"
#include "Poincare.H"
#include "Tree.H"
#include "Data_Read.H"
#include "Exception.H"
#include <iomanip>

using namespace APACIC;
using namespace ATOOLS;

Timelike_Kinematics::Timelike_Kinematics(ATOOLS::Jet_Finder *const jf,
					 double pt2min) : 
  m_pt2min(pt2min), m_pt_scheme(1), m_zrange_scheme(1), m_type(-1),
  p_jf(jf)
{
  if (rpa.gen.Beam1().IsLepton() && rpa.gen.Beam2().IsLepton())         
    m_type = 1;
  else if (!rpa.gen.Beam1().IsLepton() && !rpa.gen.Beam2().IsLepton())  
    m_type = 4;
  else {
    msg.Error()<<"ERROR in Timelike_Kinematics : "<<std::endl
	       <<"   DIS is not yet implemented in the jetfinder, "
	       <<"continue with hh-mode."<<std::endl;
    m_type = 4;
    THROW(not_implemented,"DIS is not implemented yet");
  }
}


//-----------------------------------------------------------------------
//------------------- Checks for kinematics : The shuffles --------------
//----------------------------------------------------------------------- 

bool Timelike_Kinematics::CheckZRange(Knot const * const  mo,
				      Flavour const * const d1_flav,
				      Flavour const * const d2_flav) const
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

  double z(CalcZShift(mo->z,t,t1,t2,t01,t02)), 
    e12(z*z*mo->E2), e22((1.-z)*(1.-z)*mo->E2);
  bool do1(false), do2(false);
  if (d1->stat>0 && !CheckZRange(d1->z,e12,t1,sqr(d1_flav[0].PSMass()),
				 sqr(d1_flav[1].PSMass()))) do1=true;
  if (d2->stat>0 && !CheckZRange(d2->z,e22,t2,sqr(d2_flav[0].PSMass()),
				 sqr(d2_flav[1].PSMass()))) do2=true;

  if (!do1 && !do2) return true;
  else if (!do1 && do2) d2->stat=3;
  else if (do1 && !do2) d1->stat=3;
  else if (d1->t>d2->t) d1->stat=3;
                   else d2->stat=3;
  return false;
}

bool Timelike_Kinematics::
CheckZRange(const double z, const double E2, const double t, 
	    const double t1, const double t2) const
{
  double x_p_mom(sqrt(1.-t/E2));   
  double mean_z(0.),delta_z(1.);              
  
  switch (m_zrange_scheme) {
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
CalcZShift(const double z, const double t, 
	   const double t1, const double t2, 
	   const double t01, const double t02) const
{
  double lambda1(sqr(t-t1-t2)-4.*t1*t2); 
  if (lambda1<0) {
    msg_Tracking()<<"Timelike_Kinematics::CalcZShift(..): "
		  <<"kinematics does not fit : "
		  <<t<<" -> "<<t1<<"+"<<t2<<std::endl;
    return -1.;
  }
  switch (m_zrange_scheme) {
  case 1:
    return ((2.*z-1.)*sqrt(lambda1) + (t+t1-t2))/(2.*t);
  case 0:
    double lambda0(sqr(t-t01-t02)-4.*t01*t02); 
    return  (z-(t+t01-t02)/(2.*t))*sqrt(lambda1/lambda0) + (t+t1-t2)/(2.*t);
  }
  return -1.0;
}

double Timelike_Kinematics::
CalcZShift(const double z, const double t, 
	   const double t1, const double t2) const
{
  switch (m_zrange_scheme) {
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
  bool shuffle;
  if (first) shuffle=ShuffleMoms(mo);
  else shuffle=ShuffleZ(mo);
  msg_Debugging()<<"}\n";
  return shuffle;
} 

int Timelike_Kinematics::ShuffleZ(Knot * const mo) const
{
  double t(mo->t), t1(mo->left->t), t2(mo->right->t);
  if (t - (t1+t2+2.*sqrt(t1*t2)) < rpa.gen.Accu()) {
    msg_Debugging()<<"Timelike_Kinematics::ShuffleZ(): "
		   <<"not enough virtuality: "
		   <<sqrt(t)<<" < "<<sqrt(t1)<<" + "<<sqrt(t2)<<std::endl;
    return 0; 
  }

  double t01(mo->left->tout), t02(mo->right->tout), z(mo->z);
  mo->z   = CalcZShift(z,t,t1,t2,t01,t02);
  if (FailedKinCheck(0,mo)) {
    mo->z = z;
    return 0;
  }
  mo->left->E2  = mo->z*mo->z*mo->E2;
  mo->right->E2 = (1.-mo->z)*(1.-mo->z)*mo->E2;
  msg_Debugging()<<"z = "<<mo->z<<", E_1 = "<<sqrt(mo->left->E2)
		 <<", E_2 = "<<sqrt(mo->right->E2)<<"\n";
  return 1;
}

int Timelike_Kinematics::ShuffleMoms(Knot * const mo) const
{ 
  Knot * d1(mo->left), * d2(mo->right);
  double t1(d1->t),      t2(d2->t);

  double newt(mo->part->Momentum().Abs2());
  if (dabs((mo->t - newt)/mo->t)>1.e-7) {
    msg.Error()<<"WARNING in Timelike_Kinematics::ShuffleMoms : \n"
	       <<"    Large mass deviation "<<mo->t<<" vs. "
	       <<newt<<std::endl;
  }
  double t = mo->t = newt;
  if (t - (t1+t2+2.*sqrt(t1*t2))<rpa.gen.Accu()) {
    msg_Debugging()<<"missing mass\n";
    return 0;
  }

  double r1(0.), r2(0.), z(mo->z);
  Vec4D p1(d1->part->Momentum()), p2(d2->part->Momentum());

  if ((p1.Abs2()/d1->E2 > rpa.gen.Accu()) || 
      (p2.Abs2()/d2->E2 > rpa.gen.Accu())) {
    double t1n(p1.Abs2()), t2n(p2.Abs2());
    double A(((t2 - t2n) - (t1 -t1n)) / (t + t1n - t2n));
    double B((t + t2n-t1n)/(t + t1n -t2n));
    double C(t - t1n - t2n);
    double D((2.*t2n - 2. * A*B*t1n + (A - B)*C)/
	     (2.*(t2n + B*B*t1n - B*C)));
    double E((t2n - t2 + A*A*t1n + A*C)/(t2n + B*B*t1n - B*C));
    r2 = D - sqrt(D*D- E);
    r1 = A + r2*B;
    Vec4D p1a( (1.-r1)*p1 + r2*p2 );
    Vec4D p(mo->part->Momentum());
    mo->z = p1a[0]/p[0];
  }
  else {
    double lambda(sqrt(sqr(t-t1-t2)-4.*t1*t2)); 
    r1     = (t+t2-t1-lambda)/(2.*t);
    r2     = (t-t2+t1-lambda)/(2.*t);
    mo->z  = z - r1*z + r2*(1.-z);
  } 
  if (dabs(mo->z-z) < rpa.gen.Accu()) {
    msg_Debugging()<<"shift unnecessary\n";
    // boost daughters if existent
    BoostDaughters(mo);
    // update daughter E2,z 
    Tree::UpdateDaughters(mo);
    return 1;
  }
  if (FailedKinCheck(1,mo)) {
    msg_Debugging()<<"kinematics check failed "<<mo->z<<" "<<z<<std::endl;
    mo->z = z;
    return 0;
  }
  msg_Debugging()<<"p-p_1-p_2, old = "<<(mo->part->Momentum()-
					 d1->part->Momentum()-
					 d2->part->Momentum())<<"\n";
  d1->part->SetMomentum( (1.-r1)*p1 + r2*p2 );
  d2->part->SetMomentum( (1.-r2)*p2 + r1*p1 );
  d1->E2   = mo->z*mo->z*mo->E2;
  d2->E2   = (1.-mo->z)*(1.-mo->z)*mo->E2;
  msg_Debugging()<<"           new = "<<(mo->part->Momentum()-
					 d1->part->Momentum()-
					 d2->part->Momentum())<<"\n";
  msg_Debugging()<<"z = "<<mo->z<<", p_1 = "<<d1->part->Momentum()
		 <<", p_2 = "<<d2->part->Momentum()<<"\n";
  // boost daughters if existent
  BoostDaughters(mo);
  // update daughter E2,z 
  Tree::UpdateDaughters(mo);
  return 1;
}

void Timelike_Kinematics::UpdateDaughters(Knot *const mo,
					  const bool force) const
{
//   msg_Debugging()<<METHOD<<"("<<mo->kn_no<<","<<force<<"): "
// 		 <<mo->part->Momentum()<<"\n";
//   msg_Indent();
  if (mo->left) {
    mo->left->E2=sqr(mo->z)*mo->E2;
    mo->right->E2=sqr((1.-mo->z))*mo->E2;
    DoSingleKinematics(mo,force);
    UpdateDaughters(mo->left);
    UpdateDaughters(mo->right);
  }
}


bool Timelike_Kinematics::
FailedKinCheck(const int first,Knot * const mo) const
{
  Knot * d1(mo->left), * d2(mo->right);
  if ((d1==0) || (d2==0)) return false;
  
  double E12(mo->z*mo->z*mo->E2), E22((1.-mo->z)*(1.-mo->z)*mo->E2);
  double t1(d1->t), t2(d2->t);
  // timelike daughters ?
  if (t1>E12 || t2>E22) return true;

  double p1p2(sqrt((E12-t1)*(E22-t2)));
  // triangular three momementum relation     
  if (mo->E2-mo->t - (E12-t1 + E22-t2 + 2.*p1p2) > 
      first*rpa.gen.Accu() ) return true;

  double cosreal((2.*mo->z*(1.-mo->z)*mo->E2-mo->t+t1+t2)/(2.*p1p2)); 
  // physical opening angle
  if ((dabs(cosreal) > 1.) && !(first)) return true;
  if (cosreal > 1.)  cosreal = 1.; 
  if (cosreal < -1.) cosreal = -1.; 
  mo->costh = cosreal;

  // physical deflection angle
  double coth1(-1);
  if (dabs(mo->E2-mo->t)>rpa.gen.Accu()) 
    coth1 = (mo->costh*p1p2+E12-t1)/(sqrt((mo->E2-mo->t)*(E12-t1)));
  if (dabs(coth1) > 1.+rpa.gen.Accu()) return true;
  return false;
}

//-----------------------------------------------------------------------
//--------------------- Evaluation of the kinematics --------------------
//----------------------------------------------------------------------- 
 
void Timelike_Kinematics::DoSingleKinematics(Knot * const mo,
					     const bool force) const
{
  if (!mo->left->didkin || force) {
    Vec4D p1(Vec4D(0.,0.,0.,0.)), p2(Vec4D(0.,0.,0.,0.));
    ConstructTwoVectors(mo,p1,p2);
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
  if (CheckVector(mo->part->Momentum()) || 
      CheckVector(mo->left->part->Momentum()) || 
      CheckVector(mo->right->part->Momentum())) error=1;
  if (!(mo->left->part->Momentum()+mo->right->part->Momentum()==
	mo->part->Momentum())) error=2;
  if (!(mo->left->part->Momentum()[2]<0) && 
      !(mo->left->part->Momentum()[2]>0)) error=3;
  if (error>0) {
    int op(msg.Error().precision(6));
    msg.Error()<<"Timelike_Kinematics::DoKinematics("<<mo->kn_no<<"): "
	       <<"Error "<<error<<": Momentum conservation violated.\n"
	       <<"     "<<mo->part->Momentum()<<" -> "
	       <<mo->part->Momentum().Abs2()<<" vs. "<<mo->t
	       <<" ("<<mo->part->Flav()<<","<<mo->part->Info()<<")\n"
	       <<"   = "<<mo->left->part->Momentum()<<" -> "
	       <<mo->left->part->Momentum().Abs2()<<" vs. "<<mo->left->t<<" ("
	       <<mo->left->part->Flav()<<","<<mo->left->part->Info()<<")\n"
	       <<"   + "<<mo->right->part->Momentum()<<" -> "
	       <<mo->right->part->Momentum().Abs2()
	       <<" vs. "<<mo->right->t<<" ("
	       <<mo->right->part->Flav()<<","<<mo->right->part->Info()
	       <<")\n     "<<(mo->part->Momentum()-
			      mo->left->part->Momentum()-
			      mo->right->part->Momentum())<<std::endl;
    msg.Error().precision(op);
    return false;
  }
  if (!DoKinematics(mo->left))  return false;
  if (!DoKinematics(mo->right)) return false;
  msg_Debugging()<<"p_d1 = "<<mo->left->part->Momentum()<<"\n";
  msg_Debugging()<<"p_d2 = "<<mo->right->part->Momentum()<<"\n";
  BoostDaughters(mo);
  msg_Debugging()<<"p_d1 = "<<mo->left->part->Momentum()<<"\n";
  msg_Debugging()<<"p_d2 = "<<mo->right->part->Momentum()<<"\n";
  msg_Debugging()<<"}\n";
  return true;
}


void Timelike_Kinematics::
ConstructTwoVectors(Knot *const mo,Vec4D &p1vec,Vec4D &p2vec) const
{
  double p(sqrt(mo->E2-mo->t));
  double E12(mo->left->E2), E22(mo->right->E2);
  double p1(sqrt(E12-mo->left->t)), p2(sqrt(E22-mo->right->t));
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
}

void Timelike_Kinematics::
ConstructDreiBein(Knot *const mo,Vec3D &n1,Vec3D &n2) const
{
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

bool Timelike_Kinematics::CheckVector(const Vec4D mom) const 
{
  if (mom.Abs2()>0 && mom.Abs2()<0) return 1;
  if (mom[0]<0) return 1;
  return 0;
}
 
void Timelike_Kinematics::
BoostDaughters(Vec4D pold,Vec4D pnew,const Vec4D &pmom,Knot *const mo) const
{
  msg_Debugging()<<METHOD<<"("<<(pnew-pold)<<","<<mo->kn_no<<"):\n";
  msg_Indent();
  int bigboost=0;
  Vec3D prot = cross(Vec3D(pnew),Vec3D(pold));
  Poincare bmom;  
  if (prot.Abs()/Vec3D(pnew).Abs()>100.*rpa.gen.Accu()) { 
    bigboost = 1;
    bmom = Poincare(pmom);
    Tree::BoRo(bmom,mo->left);
    Tree::BoRo(bmom,mo->right);
    bmom.Boost(pold);
    bmom.Boost(pnew);
    Vec4D dummy(pmom[0],-1.*Vec3D(pmom));
    bmom=Poincare(Vec4D(pmom[0],-1.*Vec3D(pmom)));
  }
  Vec3D na(pold);
  na=1./na.Abs() * na;
  Poincare bos1(pnew);
  Poincare bos(bos1*pold);
  Tree::BoRo(bos,mo->left);
  Tree::BoRo(bos,mo->right);
  pnew=mo->left->part->Momentum()+mo->right->part->Momentum();
  Vec3D nb(pnew);
  nb=1./nb.Abs() * nb;
  if (bigboost) {
    Tree::BoRo(bmom,mo->left);
    Tree::BoRo(bmom,mo->right);
  }
  else {
    if (!(na==nb)) {
      Poincare rot(pnew,pold);
      Tree::BoRo(rot,mo->left);
      Tree::BoRo(rot,mo->right);
    }
  }
}

void Timelike_Kinematics::BoostDaughters(Knot * const mo) const
{
  msg_Debugging()<<METHOD<<"("<<mo->kn_no<<"):\n";
  msg_Indent();
  Knot *d1(mo->left), *d2(mo->right);
  Vec4D p(d1->part->Momentum()+d2->part->Momentum());
  if (d1->left) {
    Vec4D p1old(d1->left->part->Momentum()+d1->right->part->Momentum());
    Vec4D p1new(d1->part->Momentum());
    if (p1new!=p1old && p1old!=Vec4D(0.,0.,0.,0.)) {
      BoostDaughters(p1old,p1new,p,d1);
    }
  }
  if (d2->left) {
    Vec4D p2old(d2->left->part->Momentum()+d2->right->part->Momentum());
    Vec4D p2new(d2->part->Momentum());
    if (p2new!=p2old && p2old!=Vec4D(0.,0.,0.,0.)) {
      BoostDaughters(p2old,p2new,p,d2);
    }
  }
}

//-----------------------------------------------------------------------
//------------------- Setting the colours after the shower --------------
//----------------------------------------------------------------------- 

bool Timelike_Kinematics::
ArrangeColourPartners(Particle const * const aup,
		      Knot const * const d1,Knot const * const d2) const
{
  if (!aup) return false;
  if (!d1)  return false;
  if (!d2)  return false;
  int jfmode(p_jf->Type());
  p_jf->SetType(1);
  if (p_jf->MTij2(aup->Momentum(),d1->part->Momentum()) <
      p_jf->MTij2(aup->Momentum(),d2->part->Momentum()) ) {
    p_jf->SetType(jfmode);
    return false;
  }
  p_jf->SetType(jfmode);
  return true;
}

//-----------------------------------------------------------------------
//------------------- PTs and jet checks --------------------------------
//----------------------------------------------------------------------- 

double Timelike_Kinematics::CalculateAngle(const Knot * knot) {
  return CalculateAngle(knot->t,knot->E2,knot->z);
}

double Timelike_Kinematics::
CalculateAngle(const double t,const double E2,const double z) 
{
  switch (m_angle_scheme) {
  default: return sqrt(dabs(t)/(z*(1.- z)*E2) );
  }
}

double Timelike_Kinematics::
CalcKt2(double z, double E2, double t, double t1, double t2) const
{
  switch (m_pt_scheme) {
  case 0: return m_pt_factor * z*(1.-z)*t;
  case 1: return m_pt_factor * (z*(1.-z)*t - (1.-z)*t1 + z*t2);
  default:
    double E12(z*z*E2), E22((1.-z)*(1.-z)*E2);
    double cosreal((2.*z*(1.-z)*E2-t+t1+t2)/(2.*sqrt((E12-t1)*(E22-t2))));
    return m_pt_factor * 2.*Min(E12,E22)*(1.-cosreal);
  }
}
