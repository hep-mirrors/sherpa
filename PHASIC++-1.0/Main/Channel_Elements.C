#include "Channel_Elements.H"
#include "Message.H"
#include "MathTools.H"
#include "Poincare.H"
#include "Random.H"

/* Two body decays .... */



using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

Channel_Elements PHASIC::CE;

double Channel_Elements::Isotropic2Weight(const Vec4D& p1,const Vec4D& p2)
{
  double massfactor = Channel_Basics::SqLam((p1+p2).Abs2(),p1.Abs2(),p2.Abs2());
  if (ATOOLS::IsZero(massfactor)) return 0.;  
  if (!(massfactor>0) && !(massfactor<0)) 
    ATOOLS::msg.Error()<<"Isotropic2Weight produces a nan!"<<endl;
  
  return 2./M_PI/massfactor;
}

void Channel_Elements::Isotropic2Momenta(Vec4D p,double& s1,double& s2,
					 Vec4D& p1,Vec4D& p2,
					 double ran1,double ran2)
{
  double s    = p.Abs2();
  double rs   = sqrt(dabs(s));
  Vec4D p1h;
  p1h[0]      = (s+s1-s2)/rs/2.;
  double p1m  = rs*Channel_Basics::SqLam(s,s1,s2)/2.;
  double ct   = 2.*ran1-1.;
  double st   = sqrt(1.-ct*ct);
  double phi  = 2.*M_PI*ran2;

  p1h = Vec4D(p1h[0],p1m*Vec3D(st*::sin(phi),st*cos(phi),ct));	
  Channel_Basics::Boost(0,p,p1h,p1);
  p2  = p+(-1.)*p1;
  s1 = Max(0.,p1.Abs2());
  s2 = Max(0.,p2.Abs2());
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
    ATOOLS::msg.Error()<<"Anisotropic2Weight produces a nan!"<<endl;

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

  if ((dabs(p1.Abs2()-s1)>1.e-5)) {  // explicitly not relative!
    ATOOLS::msg.Error()<<"Channel_Elements::Anisotropic2Momenta : Strong deviation in masses : "
			  <<"s1,p1: "<<s1<<";"<<p1.Abs2()<<" : "<<dabs(s1-p1.Abs2())<<endl;
  }
  if ((dabs(p2.Abs2()-s2)>1.e-5)) {  // explicitly not relative!
    ATOOLS::msg.Error()<<"Channel_Elements::Anisotropic2Momenta : Strong deviation in masses : "
			  <<"s2,p2: "<<s2<<";"<<p2.Abs2()<<" : "<<dabs(s2-p2.Abs2())<<endl;
  }
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

double Channel_Elements::MasslessPropWeight(double sexp,
					    double smin,double smax,
					    const double& s)
{
  if ((s<=smin) && (s>=smax)) return 0;

  double wt = 1./(pow(s,sexp)*Channel_Basics::PeakedWeight(0.,sexp,smin,smax,1));
  if (!(wt>0) && !(wt<0) && wt!=0) { 
    ATOOLS::msg.Error()<<"MasslessPropWeight produces a nan: "<<wt<<endl
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
  if (!(s>0) && !(s<0) && s!=0) 
    ATOOLS::msg.Error()<<"MasslessPropMomenta produced a nan !"<<endl;
  return s;
}


double Channel_Elements::ThresholdWeight(double mass,double smin,double smax,double s)
{
  if ((s<=smin) && (s>=smax)) return 0;
  double wt = s/((sqr(s)+pow(mass,4.)) * 
		 Channel_Basics::PeakedWeight(pow(mass,4.),1.,sqr(smin),sqr(smax),1)/2.);

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    ATOOLS::msg.Error()<<" In ThresholdWeight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<2.<<", "<<mass*mass<<" wt = "<<wt<<endl
			  <<"ThresholdWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::ThresholdMomenta(double mass,double smin,double smax,double ran)
{
  double s = sqrt(Channel_Basics::PeakedDist(pow(mass,4.),1.,sqr(smin),sqr(smax),1,ran));
  if (!(s>0) && !(s<0) && s!=0) ATOOLS::msg.Error()<<"ThresholdMomenta produced a nan !"<<endl;
  if ((s<smin) || (s>smax))     ATOOLS::msg.Error()<<"ThresholdMomenta out of bounds !"<<endl;
  return s;
}

double Channel_Elements::LLPropWeight(double sexp,double pole,
				      double smin,double smax,
				      const double& s)
{
  if ((s<=smin) && (s>=smax)) return 0;
  double wt = 1./(pow(pole-s,sexp)*Channel_Basics::PeakedWeight(pole,sexp,smin,smax,-1));

  if (!(wt>0) && !(wt<0) && wt!=0 ) {
    ATOOLS::msg.Error()<<" In LL_Weight : "<<smin<<" < "<<s<<" < "
			  <<smax<<" ^ "<<sexp<<", "<<pole<<" wt = "<<wt<<endl
			  <<"LLPropWeight produces a nan: "<<wt<<endl;
  }
  return wt;
}

double Channel_Elements::LLPropMomenta(double sexp,double pole,
				       double smin,double smax,
				       double ran)
{
  double s = Channel_Basics::PeakedDist(pole,sexp,smin,smax,-1,ran);
  if (!(s>0) && !(s<0) && s!=0) ATOOLS::msg.Error()<<"LLPropMomenta produced a nan !"<<endl;
  if ((s<smin) || (s>smax))     ATOOLS::msg.Error()<<"LLPropMomenta out of bounds !"<<endl;
  return s;
}

double Channel_Elements::MassivePropWeight(double mass,double width,int lim,
					   double smin,double smax,double s)
{
  double mass2 = mass*mass;
  double mw    = mass*width;
  if (lim==0) return mw/(M_PI*((s-mass2)*(s-mass2)+mw*mw));
  else {
    if ((s<smin) || (s>smax)) return 0.;
    double range  = smax-smin;
    double upper  = (smax-mass2)/mw;
    double lower  = (smin-mass2)/mw;
    double ymin   = atan(lower);
    double yrange = atan(range/(mw*(1.+lower*upper)));
    if (lower*upper<-1.) {
      if (upper>0) yrange = yrange + M_PI;
      if (upper<0) yrange = yrange - M_PI;
    }     

    double wt = mw/(yrange*((s-mass2)*(s-mass2)+mw*mw));

    if (!(wt>0) && !(wt<0) && wt!=0) {
      ATOOLS::msg.Error()<<"MassivePropWeight produces a nan!"<<endl;
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
    double range  = smax-smin;
    double upper  = (smax-mass2)/mw;
    double lower  = (smin-mass2)/mw;
    double ymin   = atan(lower);
    double yrange = atan(range/(mw*(1.+lower*upper)));
    if (lower*upper<-1.) {
      if (upper>0) yrange = yrange + M_PI;
      if (upper<0) yrange = yrange - M_PI;
    }     
    s = mass2+mw*tan(ran*yrange + ymin);
  }
  if (!(s>0) && !(s<0) && s!=0) 
    ATOOLS::msg.Error()<<"MasslessPropMomenta produced a nan !"<<endl;
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
  Vec4D p1inh,p1outh;
  p1inh[0]         = (s+s1in-s2in)/2./sabs;
  p1outh[0]        = (s+s1out-s2out)/2./sabs;
  double p1inmass  = sabs*Channel_Basics::SqLam(s,s1in,s2in)/2.; 
  double p1outmass = sabs*Channel_Basics::SqLam(s,s1out,s2out)/2.; 
  double a         = 1.;
  double aminct0   = aminct; 
  double ctmax_tmp = ctmax;
  if (aminctflag==1) {
    if (ATOOLS::IsZero(ctmax)) {
      if ( (sqr(p1outh[0])-s1in<0.) ) return 0.;
      ctmax_tmp = Channel_Basics::PseudoAngleCut(s1in,p1inh[0],s1out,p1outh[0]);
      if (ctmin < ctmax_tmp) return 0.;
    }
  }
  else {
    a= (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);
    if ( (1.>=a) && (a>=0.)) a = 1.;
    double ct = (p1inh[0]*p1outh[0]-p1in*p1out)/(p1inmass*p1outmass);
    aminct0   = a-ct;
  }
  if  ((aminct0 > ctmin+a-1.) || (aminct0 < ctmax_tmp+a-1.)) return 0.; 

  double wt = 2.*sabs/(-pow(aminct0,ctexp)*
			Channel_Basics::Hj1(ctexp,ctmin+a-1.,ctmax_tmp+a-1.)*p1outmass*M_PI);

  if (!(wt>0) && !(wt<0)) 
    ATOOLS::msg.Error()<<"TChannelWeight produces a nan!"<<endl;

  return wt;
}

int Channel_Elements::TChannelMomenta(Vec4D p1in,Vec4D p2in,Vec4D &p1out,Vec4D &p2out,  
				      double& s1out,double& s2out,double t_mass,
				      double ctexp,double ctmax,double ctmin,
				      double aminct,int aminctflag,double ran1,double ran2)
{
  /* Note : ct's maximal range : between ctmin = 0 and ctmax = 2 */

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
  
  double a         = 1.;
  double ctmax_tmp = ctmax;
  if (aminctflag==1) {
    if (ATOOLS::IsZero(ctmax)) {
      if ( (sqr(p1outh[0])-s1in<0.) ) return 0;
      ctmax_tmp = Channel_Basics::PseudoAngleCut(s1in,p1inh[0],s1out,p1outh[0]);
      if (ctmin < ctmax_tmp) return 0;
    }
  }
  else {
    a = (t_mass2-s1in-s1out+2.*p1outh[0]*p1inh[0])/(2.*p1inmass*p1outmass);
    if ( (1.>=a) && (a>=0.)) a = 1.;
  }
  aminct = Channel_Basics::Tj1(ctexp,ctmin+a-1.,ctmax_tmp+a-1.,ran1);                   
  double ct = a - aminct;
  double st;
  if (aminctflag==1) st = sqrt(aminct*(1.+ct)); 
                else st = sqrt(1.-sqr(ct));
  double phi = 2.*M_PI*ran2;
  p1outh     = Vec4D(p1outh[0],p1outmass*Vec3D(st*cos(phi),st*::sin(phi),ct)); 

  Vec4D help;
  Channel_Basics::Boost(1,pin,help,p1in);  
  
  Poincare Rot(p1inh,help);
  help = p1outh;
  Rot.Rotate(help);
  Channel_Basics::Boost(0,pin,help,p1out);

  p2out = pin+(-1.)*p1out;

  if (dabs(s1out-p1out.Abs2())>1.e-5) {
    ATOOLS::msg.Error()<<"Channel_Elements::TChannelMomenta : Strong deviation in masses : "
			  <<"s1,p1: "<<s1out<<";"<<p1out.Abs2()<<" : "<<dabs(s1out-p1out.Abs2())<<endl;
  }
  s1out = Max(0.,p1out.Abs2());
  s2out = Max(0.,p2out.Abs2());
}

double Channel_Elements::DiceYUniform(double tau, double * yrange, double * deltay, 
				      int fixed, double ran)
{
  double y0   = 0.5 * log(tau);
  if (fixed==1) return y0;
  if (fixed==2) return -y0;
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  double y    = ymin+(ymax-ymin)*ran;
  if (ATOOLS::IsZero(y)) y = 0.;
  if ((y<y0) || (y>-y0)){
    ATOOLS::msg.Error()<<"y out of bounds in DiceYUniform : "<<ymin<<" < "<<y<<" < "<<ymax<<endl;
  }
  return y;
}

double Channel_Elements::WeightYUniform(double tau, double * yrange, 
					double * deltay, int fixed,double y)
{
  if (fixed<3) return 1.;
  double y0   = 0.5 * log(tau);
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  if ((y<yrange[0]) || (y>yrange[1])) {
    ATOOLS::msg.Error()<<"y out of trivial bounds in CE.WeightYUniform : "
			  <<yrange[0]<<" < "<<y<<" < "<<yrange[1]<<endl;
  }
  return (ymax-ymin);
}

double Channel_Elements::DiceYCentral(double tau, double * yrange, double * deltay, 
				      int fixed, double ran)
{
  double y0   = 0.5 * log(tau);
  if (fixed==1) return y0;
  if (fixed==2) return -y0;
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  double y    = log(tan(ran*atan(exp(ymax))+(1.-ran)*atan(exp(ymin))));
  if (ATOOLS::IsZero(y)) y = 0.;
  if ((y<y0) || (y>-y0)){
    ATOOLS::msg.Error()<<"y out of bounds in DiceYCentral : "<<ymin<<" < "<<y<<" < "<<ymax<<endl;
  }
  return y;
}

double Channel_Elements::WeightYCentral(double tau, double * yrange, 
					double * deltay, int fixed,double y)
{
  if (fixed<3) return 1.;
  double y0   = 0.5 * log(tau);
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  if ((y<yrange[0]) || (y>yrange[1])) {
    ATOOLS::msg.Error()<<"y out of trivial bounds in CE.WeightYCentral : "
			  <<yrange[0]<<" < "<<y<<" < "<<yrange[1]<<endl;
  }
  //return (atan(exp(ymax))-atan(exp(ymin)))*2.*cosh(y);
  cout<<2.*(atan(exp(ymax))-atan(exp(ymin)))<<endl;
  return 2.*(atan(exp(ymax))-atan(exp(ymin)));
}

double Channel_Elements::DiceYForward(double tau, double * yrange, double * deltay, double yexp, 
				      int fixed, double ran)
{
  double y0   = 0.5 * log(tau);
  if (fixed==1) return y0;
  if (fixed==2) return -y0;
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  double y    = Channel_Basics::PeakedDist(ymax-deltay[0],yexp,ymin,ymax,-1,ran);
  if (ATOOLS::IsZero(y)) y = 0.;
  if ((y<y0) || (y>-y0)){ 
    ATOOLS::msg.Error()<<"y out of bounds in DiceYForward : "<<ymin<<" < "<<y<<" < "<<ymax<<endl;
  }
  return y;
}

double Channel_Elements::WeightYForward(double tau, double * yrange, double * deltay, 
					double yexp, int fixed, double y)
{
  if (fixed<3) return 1;
  double y0   = 0.5 * log(tau);
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  return Channel_Basics::PeakedWeight(ymax-deltay[0],yexp,ymin,ymax,-1) * pow(ymax-deltay[0]-y,yexp);
}

double Channel_Elements::DiceYBackward(double tau, double * yrange, double * deltay, 
				       double yexp, int fixed, double ran)
{
  double y0   = 0.5 * log(tau);
  if (fixed==1) return y0;
  if (fixed==2) return -y0;
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  double y= -Channel_Basics::PeakedDist(-ymin-deltay[1],yexp,-ymax,-ymin,-1,ran);
  if ((y<y0) || (y>-y0)){ 
    ATOOLS::msg.Error()<<"y out of bounds in DiceYBackward : "<<ymin<<" < "<<y<<" < "<<ymax<<endl;
  }
  return y;
}

double Channel_Elements::WeightYBackward(double tau, double * yrange, double * deltay, 
					 double yexp, int fixed, double y)
{
  if (fixed<3) return 1;
  double y0   = 0.5 * log(tau);
  double ymax = ATOOLS::Min(yrange[1],-y0+deltay[0]);
  double ymin = ATOOLS::Max(yrange[0],y0-deltay[1]);
  return Channel_Basics::PeakedWeight(-ymin-deltay[1],yexp,-ymax,-ymin,-1) * pow(-ymin-deltay[1]+y,yexp);
}

