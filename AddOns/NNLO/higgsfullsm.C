#include "higgsfullsm.H"
#include "Tools.H"

namespace SHNNLO {

static inline double ABS(Complex x) {return sqrt(real(x)*real(x)+imag(x)*imag(x));}
static inline double ABS(double x) {return std::abs(x);}
static inline double MAX(double x, double y) {return std::max(x,y);}

static Complex cLi2(Complex x){
 const double PISQ6 = 1.644934066848226;
 double x_0 = -0.30;
 double x_1 = 0.25;
 double x_2 = 0.51;
 if (x == Complex(1.0,0.0)) return PISQ6;
 if (real(x) >= x_2) { return PISQ6 - cLi2(1.-x) - log(x)*log(1.-x) ; }
 if ((ABS(imag(x)) > 1.) || (real(x)*real(x) + imag(x)*imag(x) > 1.2))
   return - cLi2(1./x) - 0.5 * log(-x) * log(-x) - PISQ6 ;
 if (real(x) <= x_0){
   Complex zz = log(1.-x);
   return -cLi2(-x/(1.-x)) - zz*zz/2. ; }
 else if (real(x) < x_1){
   Complex z = - log(1.-x);
   Complex temp = z*(1.-z/4.*(1.-z/9.*(1.-z*z/100.
                  *(1.-5.*z*z/294.*(1.-7.*z*z/360.
                  *(1.-5.*z*z/242.*(1.-7601.*z*z/354900.
                  *(1.-91.*z*z/4146.*(1.-3617.*z*z/161840.)
                   ))))))));
   return temp; }
 else return - cLi2(-x) + cLi2(x*x)/2. ;
}

static Complex cLi2x(Complex x, Complex cosy) {
        if (x==Complex(0.0,0.0)) return 0.;
        if (cosy==Complex(1.0,0.0)) return cLi2(x);
        if (x==Complex(1.0,0.0)) return -cLi2(Complex(1.,0.))/2.-pow(log(-cosy-sqrt(-1.+cosy*cosy)),2)/4.;
        Complex sqrty = sqrt(-1.+cosy*cosy);
        return ( log(x)*(log(1.-cosy*x+sqrty*x)+log(1.-(cosy+sqrty)*x)-log(1.-2.*cosy*x+x*x)) + cLi2(x/(cosy+sqrty)) + cLi2((cosy+sqrty)*x))/2.;
}

static Complex HjI3(double s, double t, double u, double v, double mf) {
	const Complex I=Complex(0.,1.);
	const double Pi=3.141592653589793;
        Complex bbm1=4.*mf*mf*t/(s*u);
        Complex beta=(1.+sqrt(1.+bbm1))/2.;
        Complex betam1=beta-1.;
        if (ABS(bbm1)<1e-10) betam1=(mf*mf*t)/(s*u);
        Complex aap1=4.*mf*mf/v;
        Complex gamma=(1.+sqrt(1.-aap1))/2.;
        Complex gammam1=gamma-1.;
        if (ABS(aap1)<1e-10) gammam1=-mf*mf/v; 
        if (v<0.) {
                return (2.*(-cLi2((betam1 - gammam1)/betam1) + cLi2((betam1 - gammam1)/beta) + cLi2(gammam1/(-1. + beta + gamma)) - 
               cLi2(gamma/(-1. + beta + gamma)) + (-pow(log(betam1),2) + pow(log(beta),2))/2. + 
               log(-1. + gamma)*log(betam1/(-1. + beta + gamma)) + log(gamma)*log((-1. + beta + gamma)/beta)))/(-1. + 2.*beta);
        }
        else if (v<=4.*mf*mf) {
                Complex r2=(s*u)/(s*u+t*v);
                Complex sqrtr2=sqrt(r2);
                //Complex phi=acos(((-1.+(aap1-1.)+2.*beta)*sqrtr2)/aap1);
                //Complex theta=acos(((1.+(aap1-1.)-2.*beta)*sqrtr2)/aap1);
                Complex cphi=(aap1+2.*betam1)*sqrtr2/aap1,phi;
                if (imag(cphi)==0. && real(cphi)<=1. && real(cphi)>=-1. ) phi=acos(cphi.real());
                else phi = -I*log(cphi+I*sqrt(1.-cphi*cphi));
                Complex ctheta=(aap1-2.*beta)*sqrtr2/aap1,theta;
                if (imag(ctheta)==0. && real(ctheta)<=1. && real(ctheta)>=-1. ) theta=acos(ctheta.real());
                else theta = -I*log(ctheta+I*sqrt(1.-ctheta*ctheta));
                return (2.*((phi - theta)*(phi - Pi + theta) - 2.*cLi2x(sqrtr2,cphi) + 2.*cLi2x(sqrtr2,ctheta)))/(-1. + 2.*beta);
        }
        else {
                return (2.*(-cLi2((gammam1)/(-beta + gamma)) + cLi2(gamma/(-betam1 + gammam1)) + cLi2(gammam1/(-1. + beta + gamma)) - 
               cLi2(gamma/(-1. + beta + gamma)) - I*Pi*log((-1. + beta + gamma)/(betam1 - gammam1)) + 
               log(-gamma/gammam1)*log((-1. + beta + gamma)/(betam1 - gammam1))))/(-1. + 2.*beta);
        }
}

static Complex HjW3(double s, double t, double u, double v, double mf) {
	return HjI3(s, t, u, v, mf) - HjI3(s, t, u, s, mf) - HjI3(s, t, u, u, mf);
}

static Complex HjW2(double s, double mf) {
        if (s<=0)
                return 4.*pow(asinh(sqrt(-s)/(2.*mf)),2);
        else if (s<=4.*mf*mf)
                return -4.*pow(asin(sqrt(s)/(2.*mf)),2);
        else {
                const Complex ipi=Complex(0.,3.141592653589793);
                const double pisq=9.869604401089358;
                double achsm=acosh(sqrt(s)/(2*mf));
                return 4.*achsm*achsm-pisq-4.*ipi*achsm;
        }
}

static Complex HjW1(double s, double mf) {
        if (s<=0)
                return 2.*sqrt(1.-4.*mf*mf/s)*asinh(sqrt(-s)/(2.*mf));
        else if (s<=4.*mf*mf)
                return 2.*sqrt(4.*mf*mf/s-1.)*asin(sqrt(s)/(2.*mf));
        else {
                const Complex ipi=Complex(0.,3.141592653589793);
                return sqrt(1.-4.*mf*mf/s)*(2.*acosh(sqrt(s)/(2.*mf))-ipi);
        }
}

static Complex Hjb2(double s, double t, double u, double mf, double mh) {
        double Mhsq=mh*mh;
        return mf*mf/Mhsq/Mhsq*(s*(-s+u)/(s+u)+2.*t*u*(2.*s+u)/(s+u)/(s+u)*(-HjW1(Mhsq,mf)+HjW1(t,mf))+(t*u/s*(HjW2(Mhsq,mf)-2.*HjW2(t,mf)))/2.+ 
   s*s*(2.*mf*mf/(s+u)-1./2.)/(s+u)*(-HjW2(Mhsq,mf)+HjW2(t,mf))+
   (-s/4.+mf*mf)*(HjW2(Mhsq,mf)/2.+HjW2(s,mf)/2.-HjW2(t,mf)+HjW3(s,t,u,Mhsq,mf)) + ((s-12.*mf*mf-4*t*u/s)*HjW3(t,s,u,Mhsq,mf))/8.);
}

static Complex Hjb4(double s, double t, double u, double mf, double mh) {
        double Mhsq=mh*mh;
	return mf*mf/Mhsq*(-2./3.+(-1./4.+mf*mf/Mhsq)*(-HjW2(Mhsq,mf)+HjW2(s,mf)+HjW3(s,t,u,Mhsq,mf)));
}

static Complex HjA5(double s, double t, double u, double mf, double mh) {
        if (mf/mh<=1e-15) return 0.;
        double Mhsq=mh*mh;
        return mf*mf/Mhsq*(4.+4.*s/(u+t)*(HjW1(s,mf)-HjW1(Mhsq,mf))+(1.-4.*mf*mf/(u+t))*(HjW2(s,mf)-HjW2(Mhsq,mf)));
}

static Complex HjA4(double s, double t, double u, double mf, double mh) {
        if (mf/mh<=1e-15) return 0.;
        return Hjb4(s,t,u,mf,mh) + Hjb4(u,s,t,mf,mh) + Hjb4(t,u,s,mf,mh);
}

Complex HjA2(double s, double t, double u, double mf, double mh) {
        if (mf/mh<=1e-15) return 0.;
        Complex A2a=Hjb2(s,t,u,mf,mh);
        Complex A2b=Hjb2(s,u,t,mf,mh);
        double sotu=MAX(ABS(t),ABS(u)); if (sotu!=0) sotu=ABS(s/sotu);
        Complex Aab=(A2a!=Complex(0.,0)?1.+A2b/A2a:(A2b!=Complex(0.,0)?1.+A2a/A2b:1.));
        if (ABS(Aab)<1e-6 && sotu<1e-3 ) { // if there is huge cancellation between A2a and A2b
          double mhsq=mh*mh, mfsq=mf*mf;
          if (u!=0&&(ABS(t/u)<1e-4)) return -(mfsq*(4.*u + 4.*mfsq*HjW2(u,mf) - u*HjW2(u,mf)))/(2.*mhsq*mhsq*u*u);
          if (t!=0&&(ABS(u/t)<1e-4)) return -(mfsq*(4.*t + 4.*mfsq*HjW2(t,mf) - t*HjW2(t,mf)))/(2.*mhsq*mhsq*t*t);
          return s*s*(1./mhsq/mhsq*(-36.*mfsq/(t + u)/(-4.*mfsq + t + u)*(4.*mfsq - t - u + 2.*mfsq*HjW1(t + u,mf)) + 
       1./t/t/u/u*(2.*t*(mfsq*t*(64.*t - 9.*u) + 3.*t*mfsq*mfsq + 36.*mfsq*mfsq*mfsq - (5.*t + 9.*u)*t*t) + 
          2.*u*(mfsq*u*(-9.*t + 64.*u) + 3.*u*mfsq*mfsq + 36.*mfsq*mfsq*mfsq - (9.*t + 5.*u)*u*u) + 
          6.*t*u*(2.*u*t*t + t*t*t + 2.*t*u*u - 6.*mfsq*(3.*t*u + t*t + u*u) + u*u*u)/(t + u) - 
          3.*t*(2.*mfsq*t*(14.*t - 3.*u) + 4.*t*mfsq*mfsq + 24.*mfsq*mfsq*mfsq - (2.*t + 3.*u)*t*t)*HjW1(t,mf) - 
          3.*u*(2.*mfsq*u*(-3.*t + 14.*u) + 4.*u*mfsq*mfsq + 24.*mfsq*mfsq*mfsq - (3.*t + 2.*u)*u*u)*HjW1(u,mf) - 
          72.*mfsq*t*u/(t + u)*(-(mfsq*(9.*t*u + 4.*t*t + 4.*u*u)) + pow(t + u,3))/(-4.*mfsq + t + u)*HjW1(t + u,mf) - 
          18.*mfsq*(mfsq - t)*(4.*mfsq*t - t*u + 4.*mfsq*mfsq)*HjW2(t,mf) - 18.*mfsq*(mfsq - u)*(4.*mfsq*u - t*u + 4.*mfsq*mfsq)*HjW2(u,mf) + 
          1./(t + u)*(3.*(-5.*u*t*t*t*t - 2.*t*t*t*t*t - 3.*t*t*t*u*u - 3.*t*t*u*u*u - 5.*t*u*u*u*u + 
                mfsq*(50.*u*t*t*t + 28.*t*t*t*t + 36.*t*t*u*u + 50.*t*u*u*u + 28.*u*u*u*u) - 2.*u*u*u*u*u + 24.*mfsq*mfsq*mfsq*pow(t + u,2) + 
                4.*mfsq*mfsq*pow(t + u,3))*HjW1(t + u,mf) + 2.*(14.*u*t*t*t*t + 5.*t*t*t*t*t + 9.*t*t*t*u*u + 9.*t*t*u*u*u + 14.*t*u*u*u*u - 
                mfsq*(79.*u*t*t*t + 64.*t*t*t*t + 6.*t*t*u*u + 79.*t*u*u*u + 64.*u*u*u*u) + 5.*u*u*u*u*u - 36.*mfsq*mfsq*mfsq*pow(t + u,2) - 
                3.*mfsq*mfsq*pow(t + u,3) + 9.*mfsq*(4.*(t + u)*mfsq*mfsq*mfsq - mfsq*(5.*u*t*t + 4.*t*t*t + 5.*t*u*u + 4.*u*u*u) + 
                   t*u*pow(t + u,2))*HjW2(t + u,mf))))))/36.;
        }
        return A2a+A2b;
}

static Complex HA0(double mh, double mf) {
        const Complex ipi = Complex(0.,3.141592653589793);
        if (mh==0.) return 1.;
        if (mf/mh<=1e-15) return 0.;
        double z = mh*mh/mf/mf/4.;
        Complex f = 1. - ( z<1. ? (1.-z)/z*pow(asin(sqrt(z)),2) : (z-1.)/z*pow(log(sqrt(z)+sqrt(z-1.))+ipi/2.,2) ) ;
        return 3./2./z*f;
}

//----

// The followings are all normalized by results with infinite top mass

// LO ggH amplitude-squared (one quark loop induced) with full quark mass dependence
double ggH1l(double mh, double mt, double mb, double mc) {
        Complex A0 = HA0(mh,mt)+HA0(mh,mb)+HA0(mh,mc);
        double A0sq= real(A0*conj(A0));
        return A0sq;
}

// LO ggHg amplitude-squared (one quark loop induced) with full quark mass dependence
double ggHg1l(double s, double t, double u, double mh, double mt, double mb, double mc) {
        if ((ABS(t)/MAX(ABS(s),ABS(u))<1e-10)||(ABS(u)/MAX(ABS(s),ABS(t))<1e-10)||(ABS(s)/MAX(ABS(t),ABS(u))<1e-10))
          return ggH1l(mh,mt,mb,mc);
        double Msqinfmt = (1.+(pow(s,4)+pow(t,4)+pow(u,4))/pow(mh,8))/9.;
        Complex A2stu = HjA2(s,t,u,mt,mh)+HjA2(s,t,u,mb,mh)+HjA2(s,t,u,mc,mh);
        double A2stusq=real(A2stu*conj(A2stu));
        Complex A2ust = HjA2(u,s,t,mt,mh)+HjA2(u,s,t,mb,mh)+HjA2(u,s,t,mc,mh);
        double A2ustsq=real(A2ust*conj(A2ust));
        Complex A2tus = HjA2(t,u,s,mt,mh)+HjA2(t,u,s,mb,mh)+HjA2(t,u,s,mc,mh);
        double A2tussq=real(A2tus*conj(A2tus));
        Complex A4 = HjA4(s,t,u,mt,mh)+HjA4(s,t,u,mb,mh)+HjA4(s,t,u,mc,mh);
        double A4sq=(A4*conj(A4)).real();
	return (A2stusq+A2ustsq+A2tussq+A4sq)/Msqinfmt;
}

// LO qqHq amplitude-squared (one quark loop induced) with full quark mass dependence
double qqHg1l(double s, double t, double u, double mh, double mt, double mb, double mc) {
        double Msqinfmt = 4.*(u+t)*(u+t)/9./pow(mh,4);
        Complex A5 = HjA5(s,t,u,mt,mh)+HjA5(s,t,u,mb,mh)+HjA5(s,t,u,mc,mh);
        double A5sq = real(A5*conj(A5));
        return A5sq/Msqinfmt;
}

}
