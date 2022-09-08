#include "EXTRA_XS/Two2Two/NSCS_Tools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace NSCS;
using namespace ATOOLS;

double NSCS::D(const double &e3,const double &e4,const double &lam)
{
  return e3+e4-2.*e3*e4-2.*(1.-2.*lam)*sqrt(e3*e4*(1-e3)*(1-e4));
}

double NSCS::ID(const double &e3,const double &e4,const double &eta)
{
  return (1.-(e3+e4-2.*e3*e4-2.*eta)/(2.*sqrt(e3*e4*(1-e3)*(1-e4))))/2.;
}

double NSCS::N(const double &x3,const double &x4,const double &lam)
{
  return 1.+x4*(1.-2.*x3)-2.*(1.-2.*lam)*sqrt(x4*(1-x3)*(1-x3*x4));
}

double NSCS::IN(const double &x3,const double &x4,const double &eta)
{
  return (1.-(1.+x4*(1.-2.*x3)-2.*eta)/(2.*sqrt(x4*(1-x3)*(1-x3*x4))))/2.;
}

int NSCS::ClusterNNLO(NSCS_Args &a)
{
#ifdef DEBUG__NSCS_Kinematics
  DEBUG_FUNC("Sector "<<a.m_sec);
  for (int i(0);i<4;++i)
    msg_Debugging()<<a.m_p[i]<<" "<<a.m_p[i].Abs2()<<"\n";
#endif
  double s12(2.*a.m_p[0]*a.m_p[1]), s13(2.*a.m_p[0]*a.m_p[2]);
  double s14(2.*a.m_p[0]*a.m_p[3]), s23(2.*a.m_p[1]*a.m_p[2]);
  double s24(2.*a.m_p[1]*a.m_p[3]), s34(2.*a.m_p[2]*a.m_p[3]);
  double Q2(s12+s13+s14+s23+s24+s34);
  a.m_Q=sqrt(Q2);
  a.m_x1=(s13+s23+s34)/Q2;
  a.m_x2=(s14+s24+s34)/Q2;
  double e3=a.m_e13=s13*Q2/((s12+s13+s14)*(s13+s23+s34));
  double e4=a.m_e14=s14*Q2/((s12+s13+s14)*(s14+s24+s34));
  if (a.m_sec==0) {
    if (e4<e3) a.m_sec=(e4<e3/2.)?1:2;
    else a.m_sec=(e3<e4/2.)?3:4;
  }
  if (a.m_sec==1) { a.m_x3=e3; a.m_x4=2.*e4/e3; }
  else if (a.m_sec==2) { a.m_x3=e3; a.m_x4=2.*(1.-e4/e3); }
  else if (a.m_sec==3) { a.m_x3=e4; a.m_x4=2.*e3/e4; }
  else if (a.m_sec==4) { a.m_x3=e4; a.m_x4=2.*(1.-e3/e4); }
  else THROW(fatal_error,"Invalid sector");
  a.m_x2/=a.m_x1;
  double eta34=a.m_e34=2.*s34/(sqr(a.m_x1)*a.m_x2*Q2);
  double lam(a.m_lam=ID(e3,e4,sqr(e3-e4)/eta34));
  Vec4D p3(a.m_p[2]), p4(a.m_p[3]);
  Poincare zax(a.m_p[0],Vec4D::ZVEC);
  zax.Rotate(p3);
  zax.Rotate(p4);
  a.m_phi=p3.Phi();
  double sp34(2.*sqrt(lam*(1.-lam))*(e3-e4)/D(e3,e4,lam));
  double st14(sqrt(4.*e4*(1.-e4))), cp34(sqrt(1.-sp34*sp34));
  if (2.*sqr(e3-e4)/D(e3,e4,lam)/e3>2.+2.*e4/e3*(1.-2.*e3)) cp34=-cp34;
  double cp3(cos(a.m_phi)), sp3(sin(a.m_phi));
  a.m_sphi=0;
  double ref(cp3*cp34-p4[1]/(p4[0]*st14));
  if (dabs(ref-sp3*sp34)>dabs(ref+sp3*sp34)) a.m_sphi=1;
  double pt3(p3.PPerp()), ps3(a.m_x1*a.m_Q/2.);
  a.m_b=Vec4D(0.,0.,0.,pt3/ps3)-p3[3]/ps3*p3.Perp()/pt3;
  a.m_a=Vec4D(0.,cross(Vec3D(a.m_b),Vec3D(p3))/ps3);
  zax.RotateBack(a.m_b);
  zax.RotateBack(a.m_a);
#ifdef DEBUG__NSCS_Kinematics
  msg_Debugging()<<"x_1 = "<<a.m_x1<<", x_2 = "<<a.m_x2<<"\n"
		 <<"x_3 = "<<a.m_x3<<", x_4 = "<<a.m_x4<<"\n"
		 <<"e_13 = "<<a.m_e13<<", e_14 = "<<a.m_e14<<"\n"
		 <<"\\phi = "<<a.m_phi<<", lam = "<<a.m_lam<<"\n"
		 <<"sec = "<<a.m_sec<<", sign = "<<a.m_sphi<<"\n";
#endif
  return 1;
}

int NSCS::ConstructNNLO(NSCS_Args &a,const int mode)
{
#ifdef DEBUG__NSCS_Kinematics
  DEBUG_FUNC("Sector "<<a.m_sec);
  msg_Debugging()<<"x_1 = "<<a.m_x1<<", x_2 = "<<a.m_x2<<"\n"
		 <<"x_3 = "<<a.m_x3<<", x_4 = "<<a.m_x4<<"\n"
		 <<"e_13 = "<<a.m_e13<<", e_14 = "<<a.m_e14<<"\n"
		 <<"\\phi = "<<a.m_phi<<", lam = "<<a.m_lam<<"\n"
		 <<"sec = "<<a.m_sec<<", sign = "<<a.m_sphi<<"\n";
#endif
  Vec4D p1(a.m_p[0]!=Vec4D()?a.m_p[0]:Vec4D::ZVEC);
  if (p1.PPerp2()<rpa->gen.Accu()*p1[0]) p1[1]=p1[2]=0.0;
  Poincare zax(p1,Vec4D::ZVEC);
  double Q(a.m_Q), Emax(Q/2.), lam(a.m_lam);
  double x1(a.m_x1), x2(a.m_x2);
  double E3(x1*Emax), E4(E3*x2), e3, e4, r4;
  if (a.m_sec==1) { e3=a.m_x3; r4=a.m_x4/2.; e4=a.m_x3*r4; }
  else if (a.m_sec==2) { e3=a.m_x3; r4=1.-a.m_x4/2.; e4=a.m_x3*r4; }
  else if (a.m_sec==3) { e4=a.m_x3; r4=1./(a.m_x4/2.); e3=a.m_x3/r4; }
  else if (a.m_sec==4) { e4=a.m_x3; r4=1./(1.-a.m_x4/2.); e3=a.m_x3/r4; }
  else THROW(fatal_error,"Invalid sector");
  double e34(sqr(e3-e4));
  if (e34!=0.) e34/=D(e3,e4,lam);
  double E1((Q*Q-2.*Q*(E3+E4)+4.*E3*E4*e34)/
	    (2.*Q-4.*(E3*e3+E4*e4)));
  double ct13(1.-2.*e3), ct14(1.-2.*e4);
  double st13(sqrt(1.-ct13*ct13)), st14(sqrt(1.-ct14*ct14));
  double sp34(2.*sqrt(lam*(1.-lam))*(e3-e4));
  if (sp34!=0.) sp34/=D(e3,e4,lam);
  double cp34(sqrt(1.-sp34*sp34));
  if (2.*sqr(e3-e4)/D(e3,e4,lam)/e3>2.+2.*e4/e3*(1.-2.*e3)) cp34=-cp34;
  if (a.m_sphi&1) sp34=-sp34;
  double cp3(cos(a.m_phi)), sp3(sin(a.m_phi));
  a.m_p[2]=E3*Vec4D(1.,st13*cp3,st13*sp3,ct13);
  a.m_p[3]=E4*Vec4D(1.,st14*(cp3*cp34-sp3*sp34),st14*(cp3*sp34+sp3*cp34),ct14);
  a.m_p[0]=Vec4D(E1,0.,0.,E1);
  a.m_p[1]=Vec4D(Q,0.,0.,0.)-a.m_p[0]-a.m_p[2]-a.m_p[3];
  a.m_b=Vec4D(0.,-ct13*cp3,-ct13*sp3,st13);
  a.m_a=Vec4D(0.,cross(Vec3D(a.m_b),Vec3D(a.m_p[2]))/E3);
  for (int i(0);i<4;++i) zax.RotateBack(a.m_p[i]);
  zax.RotateBack(a.m_b);
  zax.RotateBack(a.m_a);
  a.m_s14s13=x2*r4;
  a.m_s34s13=x1*x2*sqr(1.-r4)/((1.-x1*(1+x2))*(1.+r4-2.*sqrt(r4)*(1.-2.*lam)));
  if (mode) {
    if (a.m_sec==1) {
      if (e4==0.0 || e3==0.0) { a.m_p[0]+=a.m_p[3]; a.m_p[3]=Vec4D();	}
      if (e3==0.0) { a.m_p[0]+=a.m_p[2]; a.m_p[2]=Vec4D(); }
    }
    else if (a.m_sec==2) {
      if (a.m_x4==0.0 || e3==0.0) { a.m_p[2]+=a.m_p[3]; a.m_p[3]=Vec4D(); }
      if (e3==0.0) { a.m_p[0]+=a.m_p[2]; a.m_p[2]=Vec4D(); }
    }
    else if (a.m_sec==3) {
      if (a.m_x4==0.0 || e4==0.0) { a.m_p[0]+=a.m_p[2]; a.m_p[2]=Vec4D(); }
      if (e4==0.0) { a.m_p[0]+=a.m_p[3]; a.m_p[3]=Vec4D(); }
    }
    else if (a.m_sec==4) {
      if (a.m_x4==0.0 || e4==0.0) { a.m_p[3]+=a.m_p[2]; a.m_p[2]=Vec4D(); }
      if (e4==0.0) { a.m_p[0]+=a.m_p[3]; a.m_p[3]=Vec4D(); }
    }
  }
#ifdef DEBUG__NSCS_Kinematics
  for (int i(0);i<4;++i)
    msg_Debugging()<<a.m_p[i]<<" "<<a.m_p[i].Abs2()<<"\n";
#endif
  return 1;
}

double NSCS::PSWeight(const NSCS_Args &a,const int mode)
{
  double Q(a.m_Q), Emax(Q/2.), W(1.0);
  if (mode&1) {
    W*=pow(Emax,4)*sqr(a.m_x1)*a.m_x2*a.m_x3/sqr(2.*M_PI);
    if (false) {//a.m_x4) {
      if (a.m_sec==1 || a.m_sec==3) {
	double F0((1.-a.m_x4/2.)/(2.*N(a.m_x3,a.m_x4/2.,a.m_lam)));
	W*=2.*F0/sqrt(a.m_lam*(1.-a.m_lam));
      }
      if (a.m_sec==2 || a.m_sec==4) {
	double F0((a.m_x4/2.)/(4.*N(a.m_x3,1.-a.m_x4/2.,a.m_lam)));
	W*=2.*F0/sqrt(a.m_lam*(1.-a.m_lam));
      }
    }
  }
  if (mode&2) {
    double E3(a.m_x1*Emax), E4(E3*a.m_x2), e3, e4;
    if (a.m_sec==1) { e3=a.m_x3; e4=a.m_x3*a.m_x4/2.; }
    else if (a.m_sec==2) { e3=a.m_x3; e4=a.m_x3*(1.-a.m_x4/2.); }
    else if (a.m_sec==3) { e4=a.m_x3; e3=a.m_x3*a.m_x4/2.; }
    else if (a.m_sec==4) { e4=a.m_x3; e3=a.m_x3*(1.-a.m_x4/2.); }
    else THROW(fatal_error,"Invalid sector");
    double e34(sqr(e3-e4));
    if (e34!=0.) e34/=D(e3,e4,a.m_lam);
    W*=(Q*Q-2.*Q*(E3+E4)+4.*E3*E4*e34)/sqr(2.*Q-4.*(E3*e3+E4*e4));
  }
  return W;
}

double NSCS::W1314(const NSCS_Args &a,const int &mode)
{
  Vec4D p1(a.m_p[0]), p2(a.m_p[1]), p3(a.m_p[2]), p4(a.m_p[3]);
  Vec4D P(p1+p2+p3+p4);
  double Q(P.Mass()), E1(p1*P/Q), E2(p2*P/Q), E3(p3*P/Q), E4(p4*P/Q);
  double r13((p1*p3)/(E1*E3)), r14((p1*p4)/(E1*E4)), r34((p3*p4)/(E3*E4));
  double r23((p2*p3)/(E2*E3)), r24((p2*p4)/(E2*E4));
  if (mode&1) r34=0.;
  double d3(r13+r23), d4(r14+r24), d3421(r34+r23+r14), d3412(r34+r13+r24);
  return (r23*r24)/(d3*d4)*(1.+r13/d3421+r14/d3412);
}
    
double NSCS::W1324(const NSCS_Args &a,const int &mode)
{
  const Vec4D &p1(a.m_p[0]), &p2(a.m_p[1]), &p3(a.m_p[2]), &p4(a.m_p[3]);
  Vec4D P(p1+p2+p3+p4);
  double Q(P.Mass()), E1(p1*P/Q), E2(p2*P/Q), E3(p3*P/Q), E4(p4*P/Q);
  double r13((p1*p3)/(E1*E3)), r14((p1*p4)/(E1*E4)), r34((p3*p4)/(E3*E4));
  double r23((p2*p3)/(E2*E3)), r24((p2*p4)/(E2*E4));
  if (mode&1) r34=0.;
  double d3(r13+r23), d4(r14+r24), d3412(r34+r13+r24);
  return (r23*r14*r34)/(d3*d4*d3412);
}

int NSCS::Theta(const NSCS_Args &a,const int mode)
{
  if (mode==1) return a.m_e14<a.m_e13/2.;
  if (mode==2) return a.m_e14>a.m_e13/2. && a.m_e14<a.m_e13;
  if (mode==3) return a.m_e13<a.m_e14/2.;
  if (mode==4) return a.m_e13>a.m_e14/2. && a.m_e13<a.m_e14;
  THROW(fatal_error,"Invalid sector");
  return -1;
}

