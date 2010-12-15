#include "POWHEG/Showers/Kinematics_Base.H"
#include "POWHEG/Tools/Singlet.H"
#include "POWHEG/Showers/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Selectors/Jet_Finder.H"

using namespace POWHEG;
using namespace ATOOLS;
using namespace std;

const Vec3D Kinematics_Base::s_ex(Vec3D(1.,0.,0.));
const Vec3D Kinematics_Base::s_ey(Vec3D(0.,1.,0.));
const Vec3D Kinematics_Base::s_ez(Vec3D(0.,0.,1.));

double Kinematics_Base::GetS
(const double &Q2,const double &y,
 const double &mi2,const double &mj2,const double &mk2) const
{
  return y*(Q2-mk2)+(1.0-y)*(mi2+mj2);
}

double Kinematics_Base::GetZ
(const double &Q2,const double &sij,const double &y,const double &zt,
 const double &mi2,const double &mk2) const
{
  double ecm=0.5*(Q2-sij-mk2), rtlam=sqr(ecm);
  if (rtlam<sij*mk2) return sqrt(-1.0);
  rtlam=sqrt(rtlam-sij*mk2);
  double gam=ecm+Sign(Q2-sij-mk2)*rtlam;
  return ecm/rtlam*(zt-mk2/dabs(gam)*(y/(1.0-y)+mi2/ecm));
}

double Kinematics_Base::GetKT2
(const double &Q2,const double &y,const double &z,
 const double &mi2,const double &mj2,const double &mk2) const
{
  return (Q2-mi2-mj2-mk2)*y*z*(1.0-z)-sqr(1.0-z)*mi2-z*z*mj2;
}

double Kinematics_Base::ConstructLN
(const double &Q2,const double &sij,
 const double &mij2,const double &mk2,
 const Vec4D &Q,Vec4D &pk,Vec4D &l,Vec4D &n) const
{
  double po=sqr(Q2-mij2-mk2)-4.0*mij2*mk2, pn=sqr(Q2-sij-mk2)-4.0*sij*mk2;
  if ((po<0.0)^(pn<0.0)) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return 0.0;
  }
  pk=(Q2+mk2-sij)/(2.0*Q2)*Q+(pk-(Q2+mk2-mij2)/(2.0*Q2)*Q)*sqrt(pn/po);
  Vec4D pij=Q-pk;
  double gam=pij*pk+Sign(Q2-sij-mk2)*sqrt(sqr(pij*pk)-sij*mk2);
  double a13=sij/gam, a2=mk2/gam, bet=1.0/(1.0-a13*a2);
  l=bet*(pij-a13*pk);
  n=bet*(pk-a2*pij);
  return gam;
}

double Kinematics_FF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2) return -1.0;
  return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
}

int Kinematics_FF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & flj,Parton *&pc)
{
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;
  
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));

  Poincare cms(p1+p2);
  cms.Boost(rp1);
  Poincare zrot(rp1,Vec4D::ZVEC);
  if (n_perp.PSpat2()<=1.0e-6) {
    msg_Debugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(0.0,cross(Vec3D(rp1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();

  double z = split->ZTest(), y = split->YTest(), phi = split->Phi();
  double mi2=sqr(p_ms->Mass(fli)),mj2=sqr(p_ms->Mass(flj));
  double mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  double mk2 = sqr(p_ms->Mass(spect->GetFlavour()));
  Vec4D q1,q2=p2,q3,Q=p1+p2,l,n;

  double Q2=Q.Abs2();
  y=GetY(Q2,split->KtTest(),z,mi2,mj2,mk2);
  double sij=GetS(Q2,y,mi2,mj2,mk2);
  double gam=ConstructLN(Q2,sij,mij2,mk2,Q,q2,l,n);
  if (gam==0.0) return -1;
  double zt=GetZ(Q2,sij,y,z,mi2,mk2);
  double ktt=GetKT2(Q2,y,zt,mi2,mj2,mk2);
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q1 = ktt*sin(phi)*l_perp;
  cms.BoostBack(q1);
  q1 += zt*l + (mi2+ktt*ktt)/(gam*zt)*n + ktt*cos(phi)*n_perp;
  q3 = Q-q2-q1;
  
  if (!IsZero(sqr(((p1+p2)-(q1+q2+q3))[0]/(p1+p2)[0])) ||
      !IsZero(sqr(((p1+p2)-(q1+q2+q3)).Abs2()/(p1+p2).Abs2()))) {
    msg_Error()<<"Error in KinematicsFF::MakeKinematics "<< 
      " Four-Momentum violation "<<(p1+p2)-(q1+q2+q3)<<"\n"<<
      " Q initial : "<<p1+p2<<" -> "<<(p1+p2)*(p1+p2)<<"\n"<< 
      " Q final   : "<<q1+q2+q3<<" -> "<<(q1+q2+q3)*(q1+q2+q3)<<"\n";
    return -1;
  }
  
  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<METHOD<<"(): Negative energy in {\n  "
		  <<q1<<"\n  "<<q2<<"\n  "<<q3<<"\n}\n";
    return -1;
  }
  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(flj,q3,pst::FS);
  else pc->SetMomentum(q3);

  return 1;
}

double Kinematics_FI::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &mi2,const double &mj2,const double &ma2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=mi2+mj2+ma2) return -1.0;
  return 1.0/(1.0-(kt2/(z*(1.0-z))+mi2*(1.0-z)/z+mj2*z/(1.0-z))/(Q2-ma2-mi2-mj2));
}

int Kinematics_FI::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & flj,Parton *&pc)
{ 
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));
  n_perp*=1.0/n_perp.PSpat();

  Poincare cms(p1+p2);
  cms.Boost(rp1);
  Vec4D l_perp(0.0,cross(Vec3D(rp1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();
  
  Vec4D q1,q2=-p2,q3,Q=p1-p2,l,n;
  double z = split->ZTest(), y = split->YTest(), x=1.0-y;
  double phi = split->Phi(), Q2=Q.Abs2(), ma2 = p_ms->Mass2(spect->GetFlavour());
  double mi2=sqr(p_ms->Mass(fli)), mj2=sqr(p_ms->Mass(flj));
  double mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  
  y=1.0-(x=GetY(Q2,split->KtTest(),z,mi2,mj2,ma2));
  double xi=-z, yt=1.0-1.0/x;
  double sij=GetS(Q2,yt,mi2,mj2,ma2);
  double gam=ConstructLN(Q2,sij,mij2,ma2,Q,q2,l,n);
  if (gam==0.0) return -1;
  double zt=GetZ(Q2,sij,yt,xi,mij2,ma2);
  double ktt=GetKT2(Q2,yt,zt,mi2,mj2,ma2);
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q1 = ktt*sin(phi)*l_perp;
  cms.BoostBack(q1);
  q1 += zt*l + (mi2+ktt*ktt)/(gam*zt)*n + ktt*cos(phi)*n_perp;
  q3 = Q-q2-q1;
  q2 = -q2;

  split->SetYTest(1.0-x);

  if (!IsZero(sqr(((p1-p2)-(q1+q3-q2))[0])/(p1-p2).Abs2()) ||
      !IsZero(((p1-p2)-(q1+q3-q2)).Abs2()/(p1-p2).Abs2()))
    msg_Error()<<"Error in KinematicsFI::MakeKinematics "<< 
      " Four-Momentum violation "<<(p1-p2)-(q1+q3-q2)<<std::endl; 
  
  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<METHOD<<"(): Negative energy in {\n  "
		  <<q1<<"\n  "<<q2<<"\n  "<<q3<<"\n}\n";
    return -1;
  }

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(flj,q3,pst::FS);
  else pc->SetMomentum(q3);
  
  return 1;
}

double Kinematics_IF::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2>=ma2+mi2+mk2) return -1.0;
  return -z/(Q2-ma2-mi2-mk2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_IF::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fla,
 const ATOOLS::Flavour & fli,Parton *&pc)
{
  if (split->Kin()==1) {
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p2;
  Vec4D n_perp(0.0,cross(Vec3D(p2),Vec3D(p1)));
  n_perp*=1.0/n_perp.PSpat();

  Poincare cms(p1+p2);
  cms.Boost(rp1);
  Vec4D l_perp(0.0,cross(Vec3D(rp1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();

  Vec4D q1=-p1,q2,q3,Q=p2-p1,l,n;
  double kt2 = split->KtTest(), z = split->ZTest(), y = split->YTest();
  double mk2 = p_ms->Mass2(spect->GetFlavour()), phi = split->Phi();
  double mi2 = p_ms->Mass2(fli), ma2 = p_ms->Mass2(fla), Q2 = Q.Abs2();
  
  y=GetY(Q2,kt2,z,ma2,mi2,mk2);
  double xi=-y, yt=1.0-1.0/z;
  double sik=GetS(Q2,yt,mi2,mk2,ma2);
  double gam=ConstructLN(Q2,sik,mk2,ma2,Q,q1,l,n);
  if (gam==0.0) return -1;
  double zt=GetZ(Q2,sik,yt,xi,mk2,ma2);
  double ktt=GetKT2(Q2,yt,zt,mi2,mk2,ma2);
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q3 = ktt*sin(phi)*l_perp;
  cms.BoostBack(q3);
  q3 += zt*l + (mi2+ktt*ktt)/(gam*zt)*n + ktt*cos(phi)*n_perp;
  q2 = Q-q1-q3;
  q1 = -q1;

  if (!IsZero(sqr(((p2-p1)-(q2+q3-q1))[0])/(p2-p1).Abs2()) ||
      !IsZero(((p2-p1)-(q2+q3-q1)).Abs2()/(p2-p1).Abs2()))
    msg_Error()<<METHOD<<"(): Momentum not conserved. Difference is "
	       <<(p2-p1)-(q2+q3-q1)<<std::endl; 
  
  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<METHOD<<"(): Negative energy in {\n  "
		  <<q1<<"\n  "<<q2<<"\n  "<<q3<<"\n}\n";
    return -1;
  }
  
  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(fli,q3,pst::FS);
  else pc->SetMomentum(q3);

  return 1;
  }
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1(p1), rp2(p2);
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));
  n_perp*=1.0/n_perp.PSpat();

  Poincare cms(p1+p2);
  cms.Boost(p1);
  Vec4D l_perp(0.0,cross(Vec3D(p1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();
  p1=rp1;
    
  Vec4D q1,q2=p2,q3,Q=p2-p1,l,n;
  double z = split->ZTest(), y = split->YTest();
  double mk2 = p_ms->Mass2(spect->GetFlavour()), phi = split->Phi();
  double mi2 = p_ms->Mass2(fli), ma2 = p_ms->Mass2(fla);
  double mai2 = p1.Abs2(), Q2 = Q.Abs2();
  
  y=GetY(Q2,split->KtTest(),z,ma2,mi2,mk2);
  double xi=dabs((1.0-z)/(z-y)), yt=y/z, rf=dabs(z-y);
  double sai=GetS(Q2,yt,ma2,mi2,mk2);
  double gam=ConstructLN(Q2,sai,mai2,mk2,Q,q2,l,n);
  if (gam==0.0) return -1;
  double zt=GetZ(Q2,sai,yt,xi,mi2,mk2);
  double ktt=GetKT2(Q2,yt,zt,mi2,ma2,mk2);
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt)*rf;
  q3 = ktt*sin(phi)*l_perp;
  cms.BoostBack(q3);
  q3 += zt*rf*l + (mi2*rf*rf+ktt*ktt)/(gam*zt*rf)*n + ktt*cos(phi)*n_perp;
  q3[0]=sqrt(mi2*rf*rf+q3.PSpat2());
  q3 *= 1.0/rf;
  q1 = q2+q3-Q;

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"Error in  Kinematics_IF::MakeKinematics (past boost) "<<endl
	       <<" negative energy "<<q1<<"\n"
	       <<"                 "<<q2<<"\n"
	       <<"                 "<<q3<<"\n";
    return -1;
  }
  
  if (!IsEqual((rp2-rp1).Abs2(),(q2+q3-q1).Abs2(),sqrt(ATOOLS::Accu()))) {
    std::cout.precision(12);
    msg_Error()<<METHOD<<"(): Faulty kinematics "<<(rp2-rp1).Abs2()
	       <<" vs. "<<(q2+q3-q1).Abs2()<<" {\n"
	       <<"  old p_1 = "<<rp1<<"\n"
	       <<"      p_2 = "<<rp2<<"\n"
	       <<"  new q_1 = "<<q1<<"\n"
	       <<"      q_2 = "<<q2<<"\n"
	       <<"      q_3 = "<<q3<<"\n}"<<std::endl;
    return -1;
  }
  
  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(fli,q3,pst::FS);
  else pc->SetMomentum(q3);

  return 1;
}

double Kinematics_II::GetY(const double &Q2,const double &kt2,const double &z,
			   const double &ma2,const double &mi2,const double &mb2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=ma2+mi2+mb2) return -1.0;
  return z/(Q2-ma2-mb2-mi2)*((kt2+mi2)/(1.0-z)+(1.0-z)*ma2);
}

int Kinematics_II::MakeKinematics
(Parton *const split,const ATOOLS::Flavour & fli,
 const ATOOLS::Flavour & newfl,Parton *&pc)
{
  if (split->Kin()==1) {
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1=p1;
  Poincare cms(p1+p2);
  cms.Boost(rp1);
  Poincare zrot(rp1,Vec4D::ZVEC);
  Vec4D n_perp(0.0,1.0,1.0,0.0);
  zrot.RotateBack(n_perp);
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(0.0,cross(Vec3D(rp1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();
  
  double kt2 = split->KtTest(), z = split->ZTest(), y = split->YTest();
  
  double ma2 = sqr(p_ms->Mass(fli));
  double mai2 = sqr(p_ms->Mass(split->GetFlavour()));
  double mb2 = sqr(p_ms->Mass(spect->GetFlavour()));
  double mi2 = sqr(p_ms->Mass(newfl)), phi = split->Phi();
  Vec4D q1,q2,q3;
  
  Vec4D  Q=p1+p2;
  double Q2=Q.Abs2(), tt=Q2-mai2-mb2, t=Q2-ma2-mi2-mb2;
  if (tt*tt<4.*mai2*mb2 || tt<0.0 ||
      t*t<4.*ma2*mb2*z*z || t<0.0) return -1;
  double xi=z*(tt+sqrt(tt*tt-4.*mai2*mb2))/(t+sqrt(t*t-4.*ma2*mb2*z*z));
  y=GetY(Q2,kt2,z,ma2,mi2,mb2);
  double gamt=p1*p2+sqrt(sqr(p1*p2)-mai2*mb2), gam=gamt/xi;
  q1=(1.0-ma2*mb2/sqr(gam))/(1.0-mai2*mb2/sqr(xi*gam))/xi*
    (p1-mai2/(xi*gam)*p2)+p2*ma2/gam;
  q2=p2;
  double a1=ma2/gam, a2=mb2/gam, bet=1.0/(1.0-a1*a2);
  Vec4D l=bet*(q1-a1*q2), n=bet*(q2-a2*q1);
  double zt=(1.0+a1*a2)/(1.0-a1*a2)*((1.0-z)-y*(1.0+a2));
  double ktt=gam*y*(1.0+a1*a2)*zt-mi2-zt*zt*ma2;
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q3 = zt*l + (mi2+ktt*ktt)/(gam*zt)*n +
    + ktt*cos(phi)*n_perp + ktt*sin(phi)*l_perp;

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(newfl,q3,pst::FS);
  else pc->SetMomentum(q3);

  return 1;
  }
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1 = p1;

  Poincare cms(p1+p2);
  cms.Boost(rp1);
  Poincare zrot(rp1,Vec4D::ZVEC);
  Vec4D n_perp(0.0,1.0,1.0,0.0);
  zrot.RotateBack(n_perp);
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(0.0,cross(Vec3D(rp1),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();

  double z = split->ZTest(), y = split->YTest();
  
  Vec4D q1,q2=-p2,q3, Q=-p1-p2,l,n;
  double ma2 = sqr(p_ms->Mass(fli)), Q2=Q.Abs2();
  double mai2 = sqr(p_ms->Mass(split->GetFlavour()));
  double mb2 = sqr(p_ms->Mass(spect->GetFlavour()));
  double mi2 = sqr(p_ms->Mass(newfl)), phi = split->Phi();

  y=GetY(Q2,split->KtTest(),z,ma2,mi2,mb2);
  double xi=1.0-1.0/(z+y), yt=-y/z, rf=z+y;
  double sai=GetS(Q2,yt,ma2,mi2,mb2);
  double gam=ConstructLN(Q2,sai,mai2,mb2,Q,q2,l,n);
  if (gam==0.0) return -1;
  double zt=GetZ(Q2,sai,yt,xi,mi2,mb2);
  double ktt=GetKT2(Q2,yt,zt,mi2,ma2,mb2);
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt)*rf;
  q3 = ktt*sin(phi)*l_perp;
  cms.BoostBack(q3);
  q3 += zt*rf*l + (mi2*rf*rf+ktt*ktt)/(gam*zt*rf)*n + ktt*cos(phi)*n_perp;
  q3[0]=sqrt(mi2*rf*rf+q3.PSpat2());
  q3 *= 1.0/rf;
  q1=-Q+q2+q3;
  q2=-q2;

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(newfl,q3,pst::FS);
  else pc->SetMomentum(q3);

  return 1;
}

