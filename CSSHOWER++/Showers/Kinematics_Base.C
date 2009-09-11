#include "CSSHOWER++/Showers/Kinematics_Base.H"
#include "CSSHOWER++/Tools/Singlet.H"
#include "CSSHOWER++/Showers/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHASIC++/Selectors/Jet_Finder.H"

using namespace CSSHOWER;
using namespace ATOOLS;
using namespace std;

const Vec3D Kinematics_Base::s_ex(Vec3D(1.,0.,0.));
const Vec3D Kinematics_Base::s_ey(Vec3D(0.,1.,0.));
const Vec3D Kinematics_Base::s_ez(Vec3D(0.,0.,1.));

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
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));

  Poincare cms(p1+p2);
  cms.Boost(p1);
  cms.Boost(p2);
  Poincare zrot(p1,Vec4D::ZVEC);
  zrot.Rotate(p1);
  zrot.Rotate(p2);
  if (n_perp.PSpat2()>1.0e-6) {
    zrot.Rotate(n_perp);
  }
  else {
    msg_Debugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
  }
  Poincare xrot(n_perp,Vec4D::XVEC);
  
  double kt = sqrt(split->KtTest()), z = split->ZTest(), y = split->YTest(), phi = split->Phi();
  double mi2=sqr(p_ms->Mass(fli)),mj2=sqr(p_ms->Mass(flj)), mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  double mk2 = sqr(p_ms->Mass(spect->GetFlavour()));
  Vec4D q1,q2,q3,Q=p1+p2;

    //the general massive case, including massive spectator
    double Q2=Q.Abs2();
    y=GetY(Q2,kt*kt,z,mi2,mj2,mk2);
    double sij=y*(Q2-mk2)+(1.0-y)*(mi2+mj2);
    double po=sqr(Q2-mij2-mk2)-4.0*mij2*mk2, pn=sqr(Q2-sij-mk2)-4.0*sij*mk2;
    if (po<0.0 || pn<0.0) {
      msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
      return -1;
    }
    q2=(Q2+mk2-sij)/(2.0*Q2)*Q+(p2-(Q2+mk2-mij2)/(2.0*Q2)*Q)*sqrt(pn/po);
    Vec4D q13=Q-q2;
    double gam=q13*q2+sqrt(sqr(q13*q2)-sij*mk2);
    double a13=sij/gam, a2=mk2/gam, bet=1.0/(1.0-a13*a2);
    Vec4D l=bet*(q13-a13*q2), n=bet*(q2-a2*q13);
    double zt=(1.0+a13*a2)/(1.0-a13*a2)*(z-a2*(y/(1.0-y)+2.0*mi2/gam/(1.0+a13*a2)));
    double ktt=gam*y/(1.0-y)*(1.0+a13*a2)*zt*(1.0-zt)-sqr(1.0-zt)*mi2-zt*zt*mj2;
    if (ktt<0.0) {
      msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
      return -1;
    }
    ktt=sqrt(ktt);
    q1 = zt*l + (mi2+ktt*ktt)/(gam*zt)*n +
      + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
      + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 
    q3 = Q-q2-q1;

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"Error in  Kinematics_FF::MakeKinematics (before boost) "<<endl
		  <<" negative energy "<<q1<<"\n"
		  <<"                 "<<q2<<"\n"
		  <<"                 "<<q3<<"\n";
    return -1;
  }
  
  if (!IsZero(sqr(((p1+p2)-(q1+q2+q3))[0]/(p1+p2)[0])) || !IsZero(sqr(((p1+p2)-(q1+q2+q3)).Abs2()/(p1+p2).Abs2()))) {
    msg_Error()<<"Error in KinematicsFF::MakeKinematics "<< 
      " Four-Momentum violation "<<(p1+p2)-(q1+q2+q3)<<"\n"<<
      " Q initial : "<<p1+p2<<" -> "<<(p1+p2)*(p1+p2)<<"\n"<< 
      " Q final   : "<<q1+q2+q3<<" -> "<<(q1+q2+q3)*(q1+q2+q3)<<"\n";
    return -2;
  }
  
  xrot.RotateBack(q1);
  xrot.RotateBack(q2);
  xrot.RotateBack(q3);
  zrot.RotateBack(q1);
  zrot.RotateBack(q2);
  zrot.RotateBack(q3);
  cms.BoostBack(q1);
  cms.BoostBack(q2);
  cms.BoostBack(q3);

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"Error in  Kinematics_FF::MakeKinematics (after boost) "<<endl
		  <<" negative energy "<<q1<<"\n"
		  <<"                 "<<q2<<"\n"
		  <<"                 "<<q3<<"\n";
    return -1;
  }
  if (p_jf) {
    bool jet(true);
    jet&=p_jf->Qij2(q1,q3,q2,fli,flj)>=split->KtVeto();
    if (p_sud->HasKernel(spect->GetFlavour(),flj,cstp::FF) ||
	p_sud->HasKernel(flj,spect->GetFlavour(),cstp::FF))
      jet&=p_jf->Qij2(q2,q3,q1,spect->GetFlavour(),flj)>=split->KtVeto();
    if (p_sud->HasKernel(fli,spect->GetFlavour(),cstp::FF) ||
	p_sud->HasKernel(spect->GetFlavour(),fli,cstp::FF))
      jet&=p_jf->Qij2(q1,q2,q3,fli,spect->GetFlavour())>=split->KtVeto();
    if (jet) {
      msg_Debugging()<<"--- Jet veto ---\n\n";
      return 0;
    }
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
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));

  Poincare cms(p1+p2);
  cms.Boost(p1);
  cms.Boost(p2);
  Poincare zrot(p1,Vec4D::ZVEC);
  zrot.Rotate(p1);
  zrot.Rotate(p2);
  zrot.Rotate(n_perp);
  
  Poincare xrot(n_perp,Vec4D::XVEC);
  
  Vec4D q1,q2,q3,Q=p1-p2;
  double kt2 = split->KtTest(), z = split->ZTest(), y = split->YTest(), x=1.0-y;
  double phi = split->Phi(), Q2=Q.Abs2(), ma2 = p_ms->Mass2(spect->GetFlavour());
  double mi2=sqr(p_ms->Mass(fli)), mj2=sqr(p_ms->Mass(flj));
  double mij2 = sqr(p_ms->Mass(split->GetFlavour())); 
  
  x=GetY(Q2,kt2,z,mi2,mj2,ma2);
  double sij=-((1.0-x)*(Q2-ma2)-(mi2+mj2))/x;
  double po=sqr(Q2-mij2-ma2)-4.0*mij2*ma2, pn=sqr(Q2-sij-ma2)-4.0*sij*ma2;
  if (po<0.0^pn<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  q2=(Q2+ma2-sij)/(2.0*Q2)*Q+(p2-(Q2+ma2-mij2)/(2.0*Q2)*Q)*sqrt(pn/po);
  Vec4D q13=Q+q2;
  double q13q2=q13*q2, gam=q13q2+Sign(q13q2)*sqrt(sqr(q13q2)-sij*ma2);
  double bet=1.0-ma2*sij/(gam*gam);
  Vec4D l=(q13-sij/gam*q2)/bet, n=(q2-ma2/gam*q13)/bet;
  double zt=(gam*gam+ma2*sij)/(gam*gam-ma2*sij)*
    (z-ma2/gam*(1.0-x+2.0*mi2/(gam+ma2*sij/gam)));
  double ktt=(gam+ma2*sij/gam)*(1.0-x)*
    zt*(1.0-zt)-sqr(1.0-zt)*mi2-zt*zt*mj2;
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q1 = zt*l + (mi2+ktt*ktt)/(gam*zt)*n
    + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
    + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 
  q3 = Q-q1+q2;

  split->SetYTest(1.0-x);

  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"Error in  Kinematics_FI::MakeKinematics (before boost) "<<endl
	       <<" negative energy "<<q1<<"\n"
	       <<"                 "<<q2<<"\n"
	       <<"                 "<<q3<<"\n";
    return -1;
  }

  if (!IsZero(sqr(((p1-p2)-(q1+q3-q2))[0])/(p1-p2).Abs2()) || !IsZero(((p1-p2)-(q1+q3-q2)).Abs2()/(p1-p2).Abs2()))
    msg_Error()<<"Error in KinematicsFI::MakeKinematics "<< 
      " Four-Momentum violation "<<(p1-p2)-(q1+q3-q2)<<std::endl; 
  
  xrot.RotateBack(q1);
  xrot.RotateBack(q2);
  xrot.RotateBack(q3);
  zrot.RotateBack(q1);
  zrot.RotateBack(q2);
  zrot.RotateBack(q3);
  cms.BoostBack(q1);
  cms.BoostBack(q2);
  cms.BoostBack(q3);
  
  if (q1[0]<0. || q2[0]<0. || q3[0]<0.) {
    msg_Tracking()<<"Error in  Kinematics_FI::MakeKinematics (past boost) "<<endl
	       <<" negative energy "<<q1<<"\n"
	       <<"                 "<<q2<<"\n"
	       <<"                 "<<q3<<"\n";
    return -1;
  }

  if (p_jf) {
    bool jet(true);
    jet&=p_jf->Qij2(q1,q3,-q2,fli,flj)>=split->KtVeto();
    if (p_sud->HasKernel(spect->GetFlavour(),flj,cstp::IF))
      jet&=p_jf->Qij2(-q2,q3,q1,spect->GetFlavour().Bar(),flj)>=split->KtVeto();
    if (p_sud->HasKernel(spect->GetFlavour(),fli,cstp::IF))
      jet&=p_jf->Qij2(-q2,q1,q3,spect->GetFlavour().Bar(),fli)>=split->KtVeto();
    if (jet) {
      msg_Debugging()<<"--- Jet veto ---\n\n";
      return 0;
    }
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
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum(), rp1(p1), rp2(p2);
  Vec4D n_perp(0.0,cross(Vec3D(p1),Vec3D(p2)));

  Poincare cms(p1+p2);
  cms.Boost(p1);
  cms.Boost(p2);
  Poincare zrot(p1,Vec4D::ZVEC);
  zrot.Rotate(p1);
  zrot.Rotate(p2);
  zrot.Rotate(n_perp);

  Poincare xrot(n_perp,Vec4D::XVEC);
    
  Vec4D q1,q2,q3,Q=p2-p1;
  double kt2 = split->KtTest(), z = split->ZTest(), y = split->YTest();
  double mk2 = p_ms->Mass2(spect->GetFlavour()), phi = split->Phi();
  double mi2 = p_ms->Mass2(fli), ma2 = p_ms->Mass2(fla);
  double mai2 = p1.Abs2(), Q2 = Q.Abs2();
  
  //the massless & massive cases
  //fix the initial state parton momentum

  y=GetY(Q2,kt2,z,ma2,mi2,mk2);
  double tt=Q2-mai2-mk2, t=Q2-ma2-mi2-mk2;
  double sik=-((1.0-z)*(Q2-ma2)-(mi2+mk2))/z;
  double xi=z*(tt-sqrt(tt*tt-4.*mai2*mk2))/
    (t-sqrt(t*t-4.*ma2*sik*z*z));
  if (tt*tt<4.*mai2*mk2 || tt>0.0 ||
      t*t<4.*ma2*sik*z*z || t>0.0) return -1;
  double p1p2=p1*p2, gamt=p1p2+Sign(p1p2)*sqrt(sqr(p1p2)-mai2*mk2);
  if (sqr(p1p2)<mai2*mk2 || IsZero(gamt,1.0e-6)) return -1;
  double bet=1.0-mai2*mk2/(gamt*gamt), gam=gamt/xi;
  Vec4D l=(p1-mai2/gamt*p2)/bet;
  Vec4D n=(p2-mk2/gamt*p1)/bet;
  l*=(1.0-mk2/gamt)/(1.0-sik/gam);
  n*=(1.0-mai2/gamt)/(1.0-ma2/gam);
  double zt=(gam*gam+ma2*sik)/(gam*gam-ma2*sik)*
    (y-ma2/gam*(1.0-z+2.0*mi2/(gam+ma2*sik/gam)));
  double ktt=(gam+ma2*sik/gam)*(1.0-z)*
    zt*(1.0-zt)-sqr(1.0-zt)*mi2-zt*zt*mk2;
  if (ktt<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return -1;
  }
  ktt=sqrt(ktt);
  q1 = l + ma2/gam*n;
  q3 = zt*n + (mi2+ktt*ktt)/(gam*zt)*l
    + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
    + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 
  q2 = Q+q1-q3;

  xrot.RotateBack(q1);
  xrot.RotateBack(q2);
  xrot.RotateBack(q3);
  zrot.RotateBack(q1);
  zrot.RotateBack(q2);
  zrot.RotateBack(q3);
  cms.BoostBack(q1);
  cms.BoostBack(q2);
  cms.BoostBack(q3);
  
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
  
  if (p_jf) {
    bool jet(true);
    jet&=p_jf->Qij2(-q1,q3,q2,split->GetFlavour(),fli,1)>=split->KtVeto();
    if (p_sud->HasKernel(spect->GetFlavour(),fli,cstp::FI) ||
	p_sud->HasKernel(fli,spect->GetFlavour(),cstp::FI))
      jet&=p_jf->Qij2(q2,q3,-q1,spect->GetFlavour(),fli)>=split->KtVeto();
    if (p_sud->HasKernel(split->GetFlavour(),spect->GetFlavour(),cstp::IF))
      jet&=p_jf->Qij2(-q1,q2,q3,split->GetFlavour(),spect->GetFlavour(),1)>=split->KtVeto();
    if (jet) {
      msg_Debugging()<<"--- Jet veto ---\n\n";
      return 0;
    }
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
  Parton * spect = split->GetSpect();
  Vec4D p1 = split->Momentum(), p2 = spect->Momentum();
  Vec4D n_perp(0.0,1.0,1.0,0.0);
  
  Poincare cms(p1+p2);
  cms.Boost(p1);
  cms.Boost(p2);
  
  Poincare zrot;
  bool rotate = false;
  if(Vec3D(p1).Sqr()!=sqr(p1[3])) {
    rotate=true;
    zrot=Poincare(p1,Vec4D::ZVEC);
    zrot.Rotate(p1);
    zrot.Rotate(p2);
    zrot.Rotate(n_perp);
  }

  Poincare xrot(n_perp,Vec4D::XVEC);
  
  double kt2 = split->KtTest(), z = split->ZTest(), y = split->YTest();
  
  //std::cout<<" kt "<<kt<<" y "<<y<<" z "<<z<<std::endl;   

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
    + ktt*cos(phi)*Vec4D(0.0,1.0,0.0,0.0) 
    + ktt*sin(phi)*Vec4D(0.0,0.0,1.0,0.0); 

  xrot.RotateBack(q1);
  xrot.RotateBack(q2);
  xrot.RotateBack(q3);
  if (rotate) {
    zrot.RotateBack(q1);
    zrot.RotateBack(q2);
    zrot.RotateBack(q3);
  }
  cms.BoostBack(q1);
  cms.BoostBack(q2);
  cms.BoostBack(q3);

  //std::cout<<"test  : "<<y<<" vs. "<<(q1*q3)/(q1*q2)<<std::endl;  

  if (p_jf) {
    bool jet(true);
    jet&=p_jf->Qij2(-q1,q3,-q2,split->GetFlavour(),newfl,1)>=split->KtVeto();
    if (p_sud->HasKernel(spect->GetFlavour(),newfl,cstp::II))
      jet&=p_jf->Qij2(-q2,q3,-q1,spect->GetFlavour().Bar(),newfl)>=split->KtVeto();
    if (jet) {
      msg_Debugging()<<"--- Jet veto ---\n\n";
      return 0;
    }
  }

  split->SetMomentum(q1);
  spect->SetMomentum(q2);
  if (pc==NULL) pc = new Parton(newfl,q3,pst::FS);
  else pc->SetMomentum(q3);
  return 1;
}

