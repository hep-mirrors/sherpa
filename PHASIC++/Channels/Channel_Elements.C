#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Elements PHASIC::CE;

void Channel_Elements::CheckMasses(const double & s1,Vec4D & p1,const double & s2,Vec4D & p2) const
{
  if (dabs((s1-p1.Abs2())/p1[0])>1.e-8) {
    msg_Error()<<METHOD<<"(): Strong deviation in masses\n"
	       <<"s1,p1: "<<s1<<";"<<p1<<" -> "<<p1.Abs2()<<" : "
	       <<dabs(s1-p1.Abs2())<<", "
	       <<"rel = "<<dabs((s1-p1.Abs2())/p1[0])<<"."<<endl;
  }
  if (dabs((s2-p2.Abs2())/p2[0])>1.e-8) {
    msg_Error()<<METHOD<<"(): Strong deviation in masses\n"
	       <<"s2,p2: "<<s2<<";"<<p2<<" -> "<<p2.Abs2()<<" : "
	       <<dabs(s2-p2.Abs2())<<", "
	       <<"rel = "<<dabs((s2-p2.Abs2())/p2[0])<<"."<<endl;
  }
}

double Channel_Elements::Isotropic2Weight(Vec4D& p1,Vec4D& p2,
					  double& ran1,double& ran2,double ctmin,double ctmax)
{
  Vec4D p1h,p=p1+p2;

  Channel_Basics::Boost(1,p,p1h,p1);
  ran1        = (p1h[3]/p1h.PSpat()-ctmin)/(ctmax-ctmin);  
  ran2        = ::asin(p1h[1]/p1h.PPerp())/(2.*M_PI);
  if(p1h[2]<0.) ran2=.5-ran2;
  if (ran2<0.) ran2+=1.;

  double massfactor = Channel_Basics::SqLam(p.Abs2(),p1.Abs2(),p2.Abs2());
  if (ATOOLS::IsZero(massfactor)) return 0.;  
  if (!(massfactor>0) && !(massfactor<0)) 
    msg_Error()<<"Isotropic2Weight produces a nan!"<<endl;
  
  return 2./M_PI/massfactor*2.0/(ctmax-ctmin);
}

double Channel_Elements::Isotropic2Weight(const Vec4D& p1,const Vec4D& p2,double ctmin,double ctmax)
{
  double massfactor = Channel_Basics::SqLam((p1+p2).Abs2(),p1.Abs2(),p2.Abs2());
  if (ATOOLS::IsZero(massfactor)) return 0.;  
  if (!(massfactor>0) && !(massfactor<0)) 
    msg_Error()<<"Isotropic2Weight produces a nan!"<<endl;
  
  return 2./M_PI/massfactor*2.0/(ctmax-ctmin);
}

void Channel_Elements::Isotropic2Momenta(Vec4D p,double s1,double s2,
					 Vec4D& p1,Vec4D& p2,
					 double ran1,double ran2,double ctmin,double ctmax)
{
  double s    = p.Abs2();
  double rs   = sqrt(dabs(s));
  Vec4D p1h;
  p1h[0]      = (s+s1-s2)/rs/2.;
  double p1m  = rs*Channel_Basics::SqLam(s,s1,s2)/2.;
  double ct   = ctmin+(ctmax-ctmin)*ran1;
  double st   = sqrt(1.-ct*ct);
  double phi  = 2.*M_PI*ran2;

  p1h = Vec4D(p1h[0],p1m*Vec3D(st*::sin(phi),st*cos(phi),ct));	
  Channel_Basics::Boost(0,p,p1h,p1);
  p2  = p+(-1.)*p1;

  CheckMasses(s1,p1,s2,p2);
}

void Channel_Elements::Isotropic2Grid(Vec4D &p1,Vec4D &p2,
				      double& ran1,double& ran2,double ctmin,double ctmax)
{
  Vec4D p1h,p=p1+p2;
  Channel_Basics::Boost(1,p,p1h,p1);
  ran1        = (p1h[3]/p1h.PSpat()-ctmin)/(ctmax-ctmin); 
  ran2        = ::asin(p1h[1]/p1h.PPerp())/(2.*M_PI);
  if(p1h[2]<0.) ran2=.5-ran2;
  if (ran2<0.) ran2+=1.;
}

double Channel_Elements::Anisotropic2Weight(double ctexp,
					    double ctmin,double ctmax,
					    const Vec4D& p1,const Vec4D& p2)
{
  Vec4D  p      = p1+p2;
  double s      = p.Abs2();
  double s1     = p1.Abs2();
  double s2     = p2.Abs2();
  double pabs   = sqrt(dabs(s));
  Vec4D p1h;  
  p1h[0]        = (s+s1-s2)/pabs/2.;
  double p1mass = pabs*Channel_Basics::SqLam(s,s1,s2)/2.;
  double pmass  = sqrt(dabs(p[0]*p[0]-s)); 
  double a      = p[0]*p1h[0]/pmass/p1mass;

  if ((1.>=a) && (a>=0.)) a = 1.0000000001;
  double ct     = (pabs*p1[0]-p[0]*p1h[0])/pmass/p1mass;
  if ((ct<ctmin) || (ct>ctmax)) return 0.;

  double wt = 1./(M_PI*Channel_Basics::SqLam(s,s1,s2)/4.*
                    pow(a+ct,ctexp)*Channel_Basics::PeakedWeight(a,ctexp,ctmin,ctmax,1));
  if (!(wt>0) && !(wt<0)) 
    msg_Error()<<"Anisotropic2Weight produces a nan!"<<endl;

  return wt;
}

void Channel_Elements::Anisotropic2Momenta(Vec4D p,double s1,double s2,
					   double ctexp,
					   double ctmin,double ctmax,
					   Vec4D& p1,Vec4D& p2,
					   double ran1,double ran2)
{
  double s        = p.Abs2();
  double pabs     = sqrt(dabs(s));
  Vec4D p1h;
  p1h[0]          = (s+s1-s2)/pabs/2.;
  double p1mass   = pabs*Channel_Basics::SqLam(s,s1,s2)/2.;
  double pmass    = sqrt(dabs(p[0]*p[0]-s)); 
  double a        = p[0]*p1h[0]/pmass/p1mass;
  if ((1.>=a) && (a>=0.)) a = 1.0000000001;
  double   ct     = Channel_Basics::PeakedDist(a,ctexp,ctmin,ctmax,1,ran1);
  double st       = sqrt(1.-sqr(ct));
  double phi      = 2.*M_PI*ran2;
  p1h             = Vec4D(p1h[0],p1mass*Vec3D(st*::sin(phi),st*cos(phi),ct));	
  Vec4D pref,p1ref;
  pref            = Vec4D(p[0],0.,0.,pmass);

  Channel_Basics::Boost(0,pref,p1h,p1ref);
  Poincare Rot(pref,p);
  p1              = p1ref;
  Rot.Rotate(p1);

  p2 = p+(-1.)*p1;  

  CheckMasses(s1,p1,s2,p2);
}


double Channel_Elements::BremsstrahlungWeight(double ctexp,
					       double ctmin,double ctmax,
                  			       const Vec4D& q,const Vec4D& p1)
{
  Vec4D  p   = q+p1;
  double sp  = p.Abs2();
  double P   = Vec3D(p).Abs();
  double sq  = q.Abs2();
  double Q   = Vec3D(q).Abs();
  double ct  = Vec3D(p)*Vec3D(q)/(P*Q);
  if ((ct>ctmax) || (ct<ctmin)) return 0.;
  double p1m = sqrt(p1.Abs2());
  double ctkin = (2.*p[0]*q[0]-sq-sp+p1m*p1m)/(2.*P*Q);
  if ((0.<ctkin) && (ctkin<1.)) ctkin = 1.;
  double amct  = ctkin - ct;
  return 1./(-2.*M_PI*pow(amct,ctexp)*Channel_Basics::Hj1(ctexp,ctkin-ctmin,ctkin-ctmax));
}


void Channel_Elements::BremsstrahlungMomenta(Vec4D& p,const double p1mass,
					      const double Eq,const double sq,
					      const double ctmin,const double ctmax,
					      const double ctexp,
					      Vec4D &q, Vec4D &p1,
					      const double ran1,const double ran2)
{
  /* Decay p -> q + p1, q is space-like with energy Eq given from outside
     cos(pq) is constriained by ctmin and ctmax. */
  double sp    = p.Abs2();
  double P     = Vec3D(p).Abs();
  Vec4D  pnorm = Vec4D(1.,0.,0.,1.);
  double Q     = Vec3D(q).Abs();
  double ctkin = (2.*p[0]*Eq-sq-sp+p1mass*p1mass)/(2.*P*Q); 
  if ((0.<ctkin) && (ctkin<1.)) ctkin = 1.;
  double cth = ctkin-Channel_Basics::Tj1(ctexp,ctkin-ctmin,ctkin-ctmax,ran1);
  double sth = sqrt(1.-cth*cth);
  double cph = cos(2.*M_PI*ran2);
  double sph = sqrt(1.-cph*cph);
  Vec4D qref = Vec4D(Eq,Q*Vec3D(sth*cph,sth*sph,cth)); 
  double** rot;
  rot = new double*[3];
  short int i;
  for (i=0;i<3;i++) rot[i] = new double[3];
  Channel_Basics::Rotat(0,p,pnorm,rot);
  Channel_Basics::Rotat(1,q,qref,rot);
  for (i=0;i<3;i++) delete[] rot[i];
  delete[] rot;
  p1 = p+(-1.)*q;  
}

/* Propagators and other 1-dimensional Distributions */

double Channel_Elements::MasslessPropWeight(double sexp,double smin,double smax,double s)
{
  if ((s<=smin) && (s>=smax)) return 0;

  double wt = 1./(pow(s,sexp)*Channel_Basics::PeakedWeight(0.,sexp,smin,smax,1));
  if (!(wt>0) && !(wt<0) && wt!=0) { 
    msg_Error()<<"MasslessPropWeight produces a nan: "<<wt<<endl
			  <<"   smin,s,smax = "<<smin<<" < "<<s<<" < "<<smax
			  <<"   sexp = "<<sexp<<endl;
  }
  return wt;
}

double Channel_Elements::MasslessPropWeight(double sexp,double smin,double smax,
					    const double s,double &ran)
{
  if (s<smin||s>smax||smin==smax) {
    ran=-1.;
    return 0.;
  }
 
  double wt = 1./(pow(s,sexp)*Channel_Basics::PeakedWeight(0.,sexp,smin,smax,s,1,ran));
  if (!(wt>0) && !(wt<0) && wt!=0) { 
    msg_Error()<<"MasslessPropWeight produces a nan: "<<wt<<endl
			  <<"   smin,s,smax = "<<smin<<" < "<<s<<" < "<<smax
			  <<"   sexp = "<<sexp<<endl;
  }
  return wt;
}

double Channel_Elements::MasslessPropMomenta(double sexp,
					     double smin,double smax,
					     double ran)
{
  double s = Channel_Basics::PeakedDist(0.,sexp,smin,smax,1,ran);
  if (!(s>0) && !(s<0) && s!=0) {
    cout.precision(12);
    cout<<"MlPMom : "<<sexp<<" "<<smin<<" "<<smax<<" "<<s<<" "<<ran<<endl;
    msg_Error()<<"MasslessPropMomenta produced a nan !"<<endl;
  }
  return s;
}

double Channel_Elements::AntennaWeight(double amin,double amax,
				       const double a,double &ran)
{
  if (a<amin||a>amax||amin==amax) {
    ran=-1.;
    return 0.;
  }
 
  double wt = 1./(a*(1.-a)*Channel_Basics::BoundaryPeakedWeight(amin,amax,a,ran));
  if (!(wt>0) && !(wt<0) && wt!=0) { 
    msg_Error()<<"AntennaWeight produces a nan: "<<wt<<endl
		       <<"   amin,a,amax = "<<amin<<" < "<<a<<" < "<<amax<<endl;
  }
  return wt;
}

double Channel_Elements::AntennaMomenta(double amin,double amax,
					double ran)
{
  double a = Channel_Basics::BoundaryPeakedDist(amin,amax,ran);
  if (!(a>0) && !(a<0) && a!=0) 
    msg_Error()<<"AntennaMomenta produced a nan !"<<endl;
  return a;
}


double Channel_Elements::ThresholdWeight(double mass,double smin,double smax,double s)
{
  if ((s<=smin) && (s>=smax)) return 0;
  double wt = s/((sqr(s)+pow(mass,4.)) * 
		 Channel_Basics::PeakedWeight(pow(mass,4.),1.,sqr(smin),sqr(smax),1)/2.);

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    msg_Error()<<" In ThresholdWeight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<2.<<", "<<mass*mass<<" wt = "<<wt<<endl
			  <<"ThresholdWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::ThresholdMomenta(double mass,double smin,double smax,double ran)
{
  double s = sqrt(Channel_Basics::PeakedDist(pow(mass,4.),1.,sqr(smin),sqr(smax),1,ran));
  if (!(s>0) && !(s<0) && s!=0) msg_Error()<<"ThresholdMomenta produced a nan !"<<endl;
  if ((s<smin) || (s>smax))     msg_Error()<<"ThresholdMomenta out of bounds !"<<endl;
  return s;
}

double Channel_Elements::ThresholdWeight(double sexp,double mass,double smin,double smax,double s)
{
  //cout<<"Channel_Elements::ThresholdWeight "<<sexp<<" "<<mass<<" "<<smin<<" "<<smax<<" "<<s<<endl;
  if ((s<=smin) && (s>=smax)) return 0.;
  double sgmin=sqrt(sqr(smin)+sqr(sqr(mass)));
  double sgmax=sqrt(sqr(smax)+sqr(sqr(mass)));
//   double wt = s/(pow(sqr(s)+pow(mass,4.),0.5*(sexp+1.)) * 
// 		 Channel_Basics::PeakedWeight(pow(mass,4.),sexp,sqr(smin),sqr(smax),1)/2.);
  double wt = s/(pow(sqr(s)+pow(mass,4.),0.5*(sexp+1.)) * 
		 Channel_Basics::PeakedWeight(0.,sexp,sgmin,sgmax,1));

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    msg_Error()<<" In ThresholdWeight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<sexp<<", "<<mass*mass<<" wt = "<<wt<<endl
			  <<"ThresholdWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::ThresholdWeight(double sexp,double mass,double smin,double smax,double s,double &ran)
{
  //cout<<"Channel_Elements::ThresholdWeight "<<sexp<<" "<<mass<<" "<<smin<<" "<<smax<<" "<<s<<endl;
  if (s<smin||s>smax||smin==smax) {
    ran=-1.;
    return 0.;
  }
  
  double sg   =sqrt(sqr(s)+sqr(sqr(mass)));
  double sgmin=sqrt(sqr(smin)+sqr(sqr(mass)));
  double sgmax=sqrt(sqr(smax)+sqr(sqr(mass)));

  double wt = s/(pow(sg,(sexp+1.)) * 
		 Channel_Basics::PeakedWeight(0.,sexp,sgmin,sgmax,sg,1,ran));

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    msg_Error()<<" In ThresholdWeight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<sexp<<", "<<mass*mass<<" wt = "<<wt<<endl
			  <<"ThresholdWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::ThresholdMomenta(double sexp,double mass,double smin,double smax,double ran)
{
  if (smin>smax) return smax;
  double sgmin=sqrt(sqr(smin)+sqr(sqr(mass)));
  double sgmax=sqrt(sqr(smax)+sqr(sqr(mass)));
  double s = sqrt(sqr(Channel_Basics::PeakedDist(0.,sexp,sgmin,sgmax,1,ran))-sqr(sqr(mass)));
  if (!(s>0) && !(s<0) && s!=0) { msg_Error()<<"ThresholdMomenta produced a nan !"<<endl;
  cout<<"Channel_Elements::ThresholdMomenta "<<sexp<<" "<<mass<<" "<<sgmax-sgmin<<" "<<s<<" "<<ran<<endl;
  if (IsEqual(sgmin,sgmax)) s=0.5*(sgmin+sgmax);
  }
  if ((s<smin) || (s>smax)) {    msg_Error()<<"ThresholdMomenta out of bounds !"<<endl;
   cout<<"Channel_Elements::ThresholdMomenta "<<sexp<<" "<<mass<<" "<<smin<<" "<<smax<<" "<<s<<" "<<ran<<endl;
   if (s<smin) s=smin;
   else s=smax;
  }
  return s;
}

double Channel_Elements::LLPropWeight(double sexp,double pole,
				      double smin,double smax,
				      double s)
{
  if ((s<=smin) && (s>=smax)) return 0;
  double wt = 1./(pow(pole-s,sexp)*Channel_Basics::PeakedWeight(pole,sexp,smin,smax,-1));

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    msg_Error()<<" In LL_Weight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<sexp<<", "<<pole<<" wt = "<<wt<<endl
			  <<"LLPropWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::LLPropWeight(double sexp,double pole,
				      double smin,double smax,
				      double s,double& ran)
{
  if (s<smin||s>smax||smin==smax) {
    ran=-1.;
    return 0.;
  }
  double wt = 1./(pow(pole-s,sexp)*Channel_Basics::PeakedWeight(pole,sexp,smin,smax,s,-1,ran));

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    msg_Error()<<" In LL_Weight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<sexp<<", "<<pole<<" wt = "<<wt<<endl
			  <<"LLPropWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::LLPropMomenta(double sexp,double pole,
				       double smin,double smax,
				       double ran)
{
  double s;
//   cout<<"LLPropMomenta: "<<sexp<<" "<<pole<<" "<<smin<<" "<<smax<<" ";
  if (smin==smax) s=smax;
  else s = Channel_Basics::PeakedDist(pole,sexp,smin,smax,-1,ran);
//   cout<<s<<endl;
  if (!(s>0) && !(s<0) && s!=0) msg_Error()<<"LLPropMomenta produced a nan !"<<endl;
  if ((s<smin) || (s>smax))     msg_Error()<<"LLPropMomenta out of bounds !"<<endl;
  return s;
}

double Channel_Elements::MassivePropWeight(double mass,double width,int lim,
					   double smin,double smax,double s)
{
  double mass2 = mass*mass;
  double mw    = mass*width;
  if (lim==0) return mw/(M_PI*((s-mass2)*(s-mass2)+mw*mw));
  else {
    if ((s<smin) || (s>smax) || smin==smax) {
      //cout<<s<<" "<<smin<<" "<<smax<<endl;
      return 0.;
    }
    double upper  = (smax-mass2)/mw;
    double lower  = (smin-mass2)/mw;
    double wt=mw/((s-mass2)*(s-mass2)+mw*mw);
    wt/=atan(upper)-atan(lower);

    if (!(wt>0) && !(wt<0) && wt!=0) {
      msg_Error()<<"MassivePropWeight produces a nan!"<<endl;
    }
    return wt;
  }
}

double Channel_Elements::MassivePropWeight(double mass,double width,int lim,
					   double smin,double smax,double s,double &ran)
{
  double mass2 = mass*mass;
  double mw    = mass*width;
  if (lim==0) return mw/(M_PI*((s-mass2)*(s-mass2)+mw*mw));
  else {
    if ((s<smin) || (s>smax) || smin==smax) {
      ran=-1.;
     return 0.;
    }

    double ymax= atan((smin-mass2)/mw);
    double ymin= atan((smax-mass2)/mw);
    double y   = atan((s-mass2)/mw);
    ran = (y-ymin)/(ymax-ymin);

    double wt=mw/((s-mass2)*(s-mass2)+mw*mw);
    wt/=ymin-ymax;

    if (!(wt>0) && !(wt<0) && wt!=0) {
      msg_Error()<<"MassivePropWeight produces a nan!"<<endl;
    }
    return wt;
  }
}

double Channel_Elements::MassivePropMomenta(double mass,double width,int lim,
					    double smin,double smax,double ran)
{
  double mass2 = mass*mass;
  double mw    = mass*width;
  double s;
  if (lim==0) {
    s = mass2+mw*tan(M_PI*(ran-0.5));
  }
  else {
    double ymax=atan((smin-mass2)/mw);
    double ymin=atan((smax-mass2)/mw);
    s = mass2+mw*tan(ymin + ran*(ymax-ymin));
//     std::cout<<" smin/max "<<smin<<" "<<smax<<" "<<mass<<" "<<width<<" "<<ran<<" => "<<s<<std::endl;
  }
//   cout<<"MPMom :  "<<smin<<" "<<smax<<" "<<s<<" "<<ran<<" "<<mass<<" "<<width<<endl;
  if (!(s>0) && !(s<0) && s!=0) 
    msg_Error()<<"MassivePropMomenta produced a nan !"<<endl;
  return s;
}

double Channel_Elements::TChannelWeight(const Vec4D& p1in,const Vec4D& p2in,
					const Vec4D& p1out,const Vec4D& p2out,  
					double t_mass,double ctexp,
					double ctmax,double ctmin,
					double aminct,int aminctflag)
{
  double t_mass2   = t_mass*t_mass;
  Vec4D pin        = p1in+p2in;
  double s         = pin.Abs2(); 
  double sabs      = sqrt(dabs(s));
  double s1in      = p1in.Abs2();
  double s2in      = p2in.Abs2();
  double s1out     = p1out.Abs2();
  double s2out     = p2out.Abs2();
  if (s1out<1.e-8) s1out=0.;
  if (s2out<1.e-8) s2out=0.;
  Vec4D p1inh,p1outh;
  p1inh[0]         = (s+s1in-s2in)/2./sabs;
  p1outh[0]        = (s+s1out-s2out)/2./sabs;
  double p1inmass  = sabs*Channel_Basics::SqLam(s,s1in,s2in)/2.; 
  double p1outmass = sabs*Channel_Basics::SqLam(s,s1out,s2out)/2.; 
  double a = (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);
  if (a<=1.0) a = 1.0;
  double ct = (p1inh[0]*p1outh[0]-p1in*p1out)/(p1inmass*p1outmass);
  
  if  ((ct < ctmin) || (ct > ctmax)) {
    msg_Error()<<"TChannelWeight: bad momenta!!!! "<<ctmin<<" - "<<ctmax<<" ("<<ct<<")"<<endl;
    msg_Error()<<"1: "<<p1in<<endl;
    msg_Error()<<"2: "<<p2in<<endl;
    msg_Error()<<"3: "<<p1out<<endl;
    msg_Error()<<"4: "<<p2out<<endl;
    return 0.;
  } 
  aminct = a-ct;
  double wt = 2.*sabs/(-pow(aminct,ctexp)*
			Channel_Basics::Hj1(ctexp,a-ctmin,a-ctmax)*p1outmass*M_PI);

  if (!(wt>0) && !(wt<0)) 
    msg_Error()<<"TChannelWeight produces a nan!"<<endl;

  return wt;
}

double Channel_Elements::TChannelWeight(const Vec4D& p1in,const Vec4D& p2in,
					const Vec4D& p1out,const Vec4D& p2out,  
					double t_mass,double ctexp,
					double ctmax,double ctmin,
					double aminct,int aminctflag,
					double &ran1,double &ran2)
{
  // Note : ct's maximal range : between ctmin = -1 and ctmax = 1 
  double t_mass2   = t_mass*t_mass;
  Vec4D pin        = p1in+p2in;
  double s         = pin.Abs2(); 
  double sabs      = sqrt(dabs(s));
  double s1in      = p1in.Abs2();
  double s2in      = p2in.Abs2();
  double s1out     = p1out.Abs2();
  double s2out     = p2out.Abs2();
  if (s1out<1.e-8) s1out=0.;
  if (s2out<1.e-8) s2out=0.;
  Vec4D p1inh,p1outh;
  p1inh[0]         = (s+s1in-s2in)/2./sabs;
  double p1inmass  = sabs*Channel_Basics::SqLam(s,s1in,s2in)/2.; 
  p1inh            = Vec4D(p1inh[0],0.,0.,p1inmass);
  p1outh[0]        = (s+s1out-s2out)/2./sabs;
  double p1outmass = sabs*Channel_Basics::SqLam(s,s1out,s2out)/2.; 
  
  double a = (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);
  if (a<=1.0) a = 1.0;

  Vec4D help=p1out;
  Channel_Basics::Boost(1,pin,p1outh,help);
  help=p1in;
  Channel_Basics::Boost(1,pin,p1inh,help);  
//     if(!IsEqual(sqrt(s),pin[0])){
//       cout<<"2 bp1out="<<p1out<<"->"<<p1outh<<endl;
//       cout<<"2 bp1in= "<<p1in<<"->"<<p1inh<<endl;
//     }
  Poincare Rot(Vec4D(1.,0.,0.,1.),p1inh);
  Rot.RotateBack(p1outh);
//    cout<<" p1outh="<<p1outh<<endl;
//    cout<<" sphi/ct="<<p1outh[2]/p1outh.PPerp()<<"/"<<p1outh[3]/p1outh.PSpat()<<endl;
  
  double pa1;
  if (dabs(a-ctmax)<1.e-14) pa1 = 0.;
  else pa1 = pow(a-ctmax,1.-ctexp);
  double ct=p1outh[3]/p1outh.PSpat();
  if (ct<ctmin||ct>ctmax) {
    ran1=ran2=-1.;
    msg_Error()<<"TChannelWeight: bad momenta!!!! "<<ctmin<<" - "<<ctmax<<" ("<<ct<<")"<<endl;
    msg_Error()<<"1: "<<p1in<<endl;
    msg_Error()<<"2: "<<p2in<<endl;
    msg_Error()<<"3: "<<p1out<<endl;
    msg_Error()<<"4: "<<p2out<<endl;
    return 0.;
  }
  ran1        = (pow(a-ct,1.-ctexp)-pa1);
  ran1/=(pow(a-ctmin,1.-ctexp)-pa1); 
  ran2        = ::asin(p1outh[2]/p1outh.PPerp())/(2.*M_PI);
  if(p1outh[1]<0.) ran2=.5-ran2;
  if (ran2<0.) ran2+=1.;

  aminct = a-ct;
  double wt = 2.*sabs/(-pow(aminct,ctexp)*
			Channel_Basics::Hj1(ctexp,a-ctmin,a-ctmax)*p1outmass*M_PI);

  if (!(wt>0) && !(wt<0)) 
    msg_Error()<<"TChannelWeight produces a nan!"<<endl;

  return wt;
}

int Channel_Elements::TChannelMomenta(Vec4D p1in,Vec4D p2in,Vec4D &p1out,Vec4D &p2out,  
				      double s1out,double s2out,double t_mass,
				      double ctexp,double ctmax,double ctmin,
				      double aminct,int aminctflag,double ran1,double ran2)
{
  // Note : ct's maximal range : between ctmin = -1 and ctmax = 1 
  double t_mass2   = t_mass*t_mass;
  Vec4D pin        = p1in+p2in;
  double s         = pin.Abs2(); 
  double sabs      = sqrt(dabs(s));
  double s1in      = p1in.Abs2();
  double s2in      = p2in.Abs2();
  Vec4D p1inh,p1outh;
  p1inh[0]         = (s+s1in-s2in)/2./sabs;
  double p1inmass  = sabs*Channel_Basics::SqLam(s,s1in,s2in)/2.; 
  p1inh            = Vec4D(p1inh[0],0.,0.,p1inmass);
  p1outh[0]        = (s+s1out-s2out)/2./sabs;
  double p1outmass = sabs*Channel_Basics::SqLam(s,s1out,s2out)/2.; 
  
  double a = (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);
  if (a<=1.0) a = 1.0;
//      cout<<"TChannelMomenta"<<endl;
//      cout<<" a="<<a<<" "<<a-ctmin<<" "<<a-ctmax<<endl;
  if (dabs(a-ctmax)<1.e-14) a=ctmax;
  aminct = Channel_Basics::Tj1(ctexp,a-ctmin,a-ctmax,ran1);                   
  double ct = a - aminct;
  double st;
  if (aminctflag==1) st = sqrt(aminct*(1.+ct)); 
                else st = sqrt(1.-sqr(ct));
  double phi = 2.*M_PI*ran2;
  p1outh     = Vec4D(p1outh[0],p1outmass*Vec3D(st*cos(phi),st*::sin(phi),ct)); 
//     if(!IsEqual(sqrt(s),pin[0])){
//       cout.precision(12);
//      cout<<"1 p1outh="<<p1outh<<endl;
//      cout<<"1 sphi/ct/a="<<::sin(phi)<<"/"<<ct<<"/"<<endl;
//      cout<<"1 rans "<<ran1<<" "<<ran2<<endl;    
//     }
  Vec4D help;
  Channel_Basics::Boost(1,pin,help,p1in);  
//     if(!IsEqual(sqrt(s),pin[0])){
//       cout<<"1 bp1in= "<<p1in<<"<-"<<help<<endl;
//       cout<<"1 p1inh= "<<p1inh<<endl;
//     }
  
  Poincare Rot(p1inh,help);
  help = p1outh;
  Rot.Rotate(help);
  Channel_Basics::Boost(0,pin,help,p1out);
//     if(!IsEqual(sqrt(s),pin[0])){
//       cout<<"1 bp1out="<<p1out<<"<-"<<help<<"<-"<<p1outh<<endl;
//     }

  p2out = pin+(-1.)*p1out;

  CheckMasses(s1out,p1out,s2out,p2out);
  return 0;
}

void Channel_Elements::TChannelGrid(const Vec4D& p1in,const Vec4D& p2in,
					const Vec4D& p1out,const Vec4D& p2out,  
					double t_mass,double ctexp,
					double ctmax,double ctmin,
					double &ran1,double &ran2)
{
  // Note : ct's maximal range : between ctmin = -1 and ctmax = 1 
  double t_mass2   = t_mass*t_mass;
  Vec4D pin        = p1in+p2in;
  double s         = pin.Abs2(); 
  double sabs      = sqrt(dabs(s));
  double s1in      = p1in.Abs2();
  double s2in      = p2in.Abs2();
  double s1out     = p1out.Abs2();
  double s2out     = p2out.Abs2();
  if (s1out<1.e-8) s1out=0.;
  if (s2out<1.e-8) s2out=0.;
  Vec4D p1inh,p1outh;
  p1inh[0]         = (s+s1in-s2in)/2./sabs;
  double p1inmass  = sabs*Channel_Basics::SqLam(s,s1in,s2in)/2.; 
  p1inh            = Vec4D(p1inh[0],0.,0.,p1inmass);
  p1outh[0]        = (s+s1out-s2out)/2./sabs;
  double p1outmass = sabs*Channel_Basics::SqLam(s,s1out,s2out)/2.; 
  
  double a = (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);

  Vec4D help=p1out;
  Channel_Basics::Boost(1,pin,p1outh,help);
  help=p1in;
  Channel_Basics::Boost(1,pin,p1inh,help);  
//     if(!IsEqual(sqrt(s),pin[0])){
//       cout<<"2 bp1out="<<p1out<<"->"<<p1outh<<endl;
//       cout<<"2 bp1in= "<<p1in<<"->"<<p1inh<<endl;
//     }
  Poincare Rot(Vec4D(1.,0.,0.,1.),p1inh);
  Rot.RotateBack(p1outh);
//    cout<<" p1outh="<<p1outh<<endl;
//    cout<<" sphi/ct="<<p1outh[2]/p1outh.PPerp()<<"/"<<p1outh[3]/p1outh.PSpat()<<endl;
  
  ctexp=1.-ctexp;
  double pa1;
  if (dabs(a-ctmax)<1.e-14) pa1 = 0.;
  else pa1 = pow(a-ctmax,ctexp);
  double ct=p1outh[3]/p1outh.PSpat();
  if (ct<ctmin||ct>ctmax) {
    ran1=ran2=-1.;
    return;
  }
  ran1        = (pow(a-ct,ctexp)-pa1);
  ran1/=(pow(a-ctmin,ctexp)-pa1); 
  ran2        = ::asin(p1outh[2]/p1outh.PPerp())/(2.*M_PI);
  if(p1outh[1]<0.) ran2=.5-ran2;
  if (ran2<0.) ran2+=1.;
//     if(!IsEqual(sqrt(s),pin[0])){
//      cout<<"2 p1outh="<<p1outh<<endl;
//      cout<<"2 sphi/ct/a="<<p1outh[2]/p1outh.PPerp()<<"/"<<p1outh[3]/p1outh.PSpat()<<"/"<<a<<endl;
//      cout<<"2 rans "<<ran1<<" "<<ran2<<endl;    
//     }
}

// treated from here

double Channel_Elements::DiceYUniform(const double tau,const Double_Container &xinfo,
				      const Double_Container &yinfo,const double ran,const int mode) const
{
  /*!
    The boundaries for y are 
    \begin{align}
    \frac{1}{2}\log\frac{x_{1, min}^2}{\tau} \le y \le \frac{1}{2}\frac{x_{1, max}^2}{\tau}
    \frac{1}{2}\log\frac{\tau}{x_{2, max}^2} \le y \le \frac{1}{2}\frac{\tau}{x_{2, min}^2}
    \end{align}
    where $x_{1/2, max}$ stem from the corresponding Base or the hard process respectively and
    x_{1, min} = xinfo[0] x_{1, max} = xinfo[1]
    x_{2, min} = xinfo[2] x_{2, max} = xinfo[3]
  */
  double logtau=0.5*log(tau);
  if (mode==1) return logtau;
  if (mode==2) return -logtau;
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  double y=ymin+(ymax-ymin)*ran;
  if (ATOOLS::IsZero(y)) y=0.;
  if (y<ymin || y>ymax){
    msg_Error()<<"Channel_Elements::DiceYUniform("<<tau<<","<<xinfo<<","
		       <<yinfo<<"): "<<" Y out of bounds ! "<<std::endl<<"   ymin, ymax vs. y : "
		       <<ymin<<" "<<ymax<<" vs. "<<y<<endl;
  // If y is close to any bound, set it to this bound
    if (ATOOLS::IsEqual(y, ymin)) 
       { msg_Error()<<"Setting y to lower bound  ymin="<<ymin<<endl;
	 y = ymin; }
    if (ATOOLS::IsEqual(y, ymax)) 
      { msg_Error()<<"Setting y to upper bound ymax="<<ymax<<endl;
	 y = ymax; }
  }
  return y;
}

double Channel_Elements::WeightYUniform(const double tau,const Double_Container &xinfo,
					const Double_Container &yinfo,double& ran,const int mode) const
{
  /*
    See DiceYUniform for details
  */
  if (mode!=3) return 1.;
  double logtau=0.5*log(tau);
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  if (yinfo[2]<ymin || yinfo[2]>ymax) return 0.0;
  ran = (yinfo[2]-ymin)/(ymax-ymin);
  return (ymax-ymin);
}

const double pre=1.0;

double Channel_Elements::DiceYCentral(const double tau,const Double_Container &xinfo,
				      const Double_Container &yinfo,const double ran,const int mode) const
{
  double logtau=0.5*log(tau);
  if (mode==1) return logtau;
  if (mode==2) return -logtau;
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  double y=pre*tan(ran*atan(ymax/pre)+(1.-ran)*atan(ymin/pre));
  if (ATOOLS::IsZero(y)) y=0.;
  if (y<ymin || y>ymax){
    msg_Error()<<"Channel_Elements::DiceYCentral("<<tau<<","<<xinfo<<","
		       <<yinfo<<"): "<<" Y out of bounds ! "<<std::endl<<"   ymin, ymax vs. y : "
		       <<ymin<<" "<<ymax<<" vs. "<<y<<endl;
    if (ATOOLS::IsEqual(y, ymin)) 
       { msg_Error()<<"Setting y to lower bound  ymin="<<ymin<<endl;
	 y = ymin; }
    if (ATOOLS::IsEqual(y, ymax)) 
      { msg_Error()<<"Setting y to upper bound ymax="<<ymax<<endl;
	 y = ymax; }
  }
  return y;
}

double Channel_Elements::WeightYCentral(const double tau,const Double_Container &xinfo,
					const Double_Container &yinfo,double& ran,const int mode) const
{
  if (mode!=3) return 1.;
  double logtau=0.5*log(tau);
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  if (yinfo[2]<ymin || yinfo[2]>ymax) return 0.0;
  double atey = atan(ymin/pre);
  double wt = atan(ymax/pre)-atey;
  ran = (atan(yinfo[2]/pre)-atey)/wt;
  return wt/pre*(pre*pre+yinfo[2]*yinfo[2]);
}

double Channel_Elements::DiceYForward(const double yexponent,const double tau,
				      const Double_Container &xinfo,const Double_Container &yinfo, 
				      const double ran,const int mode) const
{
  double logtau=0.5*log(tau);
  if (mode==1) return logtau;
  if (mode==2) return -logtau;
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  double ypeak = ymax-xinfo[3];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak,ymax)) ypeak*=1.00000001;

  double y=Channel_Basics::PeakedDist(ypeak,yexponent,ymin,ymax,-1,ran);
  if (ATOOLS::IsZero(y)) y=0.;
  if (y<ymin || y>ymax){ 
    msg_Error()<<"Channel_Elements::DiceYForward("<<tau<<","<<xinfo<<","
		       <<yinfo<<"): "<<" Y out of bounds ! "<<std::endl<<"   ymin, ymax vs. y : "
		       <<ymin<<" "<<ymax<<" vs. "<<y<<endl;
     if (ATOOLS::IsEqual(y, ymin)) 
       { msg_Error()<<"Setting y to lower bound  ymin="<<ymin<<endl;
	 y = ymin; }
    if (ATOOLS::IsEqual(y, ymax)) 
      { msg_Error()<<"Setting y to upper bound ymax="<<ymax<<endl;
	 y = ymax; }
  }
  //std::cout<<ymin<<" "<<ymax<<" vs. "<<y<<endl;
  return y;
}

double Channel_Elements::WeightYForward(const double yexponent,const double tau,
					const Double_Container &xinfo,
					const Double_Container &yinfo,double& ran,const int mode) const
{
  if (mode!=3) return 1.;
  double logtau=0.5*log(tau);
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  if (yinfo[2]<ymin || yinfo[2]>ymax) return 0.0;
  double ypeak = ymax-xinfo[3];
  if (yexponent>=1. && ATOOLS::IsEqual(ypeak,ymax)) ypeak*=1.00000001;

  double wt = Channel_Basics::PeakedWeight(ypeak,yexponent,ymin,ymax,yinfo[2],-1,ran)* 
    pow(ypeak-yinfo[2],yexponent);
    if (!(wt>0) && !(wt<0) && wt!=0) {
      msg_Error()<<"WeightYForward produces a nan!"<<endl
			 <<ymax<<" "<<ymin<<" "<<yexponent<<" "<<yinfo[2]<<" "<<xinfo[3]<<endl;
      abort();
    }
  return wt;
}

double Channel_Elements::DiceYBackward(const double yexponent,const double tau,
				       const Double_Container &xinfo,const Double_Container &yinfo, 
				       const double ran,const int mode) const
{
  double logtau=0.5*log(tau);
  if (mode==1) return logtau;
  if (mode==2) return -logtau;
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  double y=-Channel_Basics::PeakedDist(-ymin-xinfo[1],yexponent,-ymax,-ymin,-1,ran);
  if (ATOOLS::IsZero(y)) y=0.;
  if (y<ymin || y>ymax){ 
    std::cout.precision(14);
    msg_Error()<<"Channel_Elements::DiceYBackward("<<tau<<","<<xinfo<<","
		       <<yinfo<<"): ";
    std::cout.precision(14);
msg_Error()<<" Y out of bounds ! "<<std::endl<<"   ymin, ymax vs. y : "
		   <<ymin<<" "<<ymax<<" vs. "<<y<<endl;
    if (ATOOLS::IsEqual(y, ymin)) 
       { msg_Error()<<"Setting y to lower bound  ymin="<<ymin<<endl;
	 y = ymin; }
    if (ATOOLS::IsEqual(y, ymax)) 
      { msg_Error()<<"Setting y to upper bound ymax="<<ymax<<endl;
	 y = ymax; }

  }
  return y;
}

double Channel_Elements::WeightYBackward(const double yexponent,const double tau,
					 const Double_Container &xinfo,
					 const Double_Container &yinfo,double& ran,const int mode) const
{
  if (mode!=3) return 1.;
  double logtau=0.5*log(tau);
  double ymin=ATOOLS::Max(xinfo[0]-logtau,logtau-xinfo[3]);
  double ymax=ATOOLS::Min(xinfo[1]-logtau,logtau-xinfo[2]);
  ymin=ATOOLS::Max(yinfo[0],ymin);
  ymax=ATOOLS::Min(yinfo[1],ymax);
  if (yinfo[2]<ymin || yinfo[2]>ymax) return 0.0;
  double wt = Channel_Basics::PeakedWeight(-ymin-xinfo[1],yexponent,-ymax,-ymin,-yinfo[2],-1,ran)* 
    pow(-ymin-xinfo[1]+yinfo[2],yexponent);
    if (!(wt>0) && !(wt<0) && wt!=0) {
      msg_Error()<<"WeightYForward produces a nan!"<<endl
			 <<ymax<<" "<<ymin<<" "<<yexponent<<" "<<yinfo[2]<<" "<<xinfo[3]<<endl;
      abort();
    }
  return wt;
}
