#include "MathTools.H"
#include <iostream>

namespace ATOOLS {

  // calculates the logarithm of the Gammafunction
  double Gammln(double xx)
  {
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    short int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
  }
  // (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. 


  double ReIncompleteGamma0(double x, double prec) {
    const double euler_gamma=GAMMA_E;

    /*
     Fuer die unvollstaendige Gammafunktion Gamma(0,x) (x reel) existieren
     Entwicklungen sowohl fuer positive als auch fuer negative x.
     Fuer positive x kann man um grosse x entwickeln, fuer negative x nur um
     kleine x:
       x>0  x->infinity
     Gamma(0,x) = - Exp[-x] * ( Sum_n=1^Infty   (n-1)! / (-x)^n)
       x->0
     Re[Gamma(0,x)] = - EulerGamma - Log[Abs[x]] - Sum_n=1^Infty(-x)^n/(n*n!)
     Im = - I Pi  fuer x<0
     Im = 0       fuer x>0
 
       z.B.
     Gamma(0,-.1) = 1.6228128139692766750 - 3.1415926535897932385 I
     Gamma(0,.1)  = 1.8229239584193906661


     Fuer positive argumente ist die funktion alternierend
   */

    double sum= -euler_gamma -log(dabs(x));
    double i  = 1;
    double ai = -x;  
    for (;;) {
      sum-=ai;
      ai*=-x*i/sqr(i+1);
      i+=1.;
      if (dabs(ai/sum)<prec) break;
      if (i>2000) {
	std::cerr<<" ERROR in ReIncompletGamma0("<<x<<")"<<std::endl;
	std::cerr<<"       "<<i<<" iteration and error="<<dabs(ai/sum)<<std::endl;
	std::cerr<<"       still bigger than wanted "<<prec<<std::endl;
	std::cerr<<"       returning "<<sum-ai<<std::endl;
      }
    }
    sum-=ai;

    return sum;

  }

  double polevl(double x,double* coef,int N )
  {
    double ans;
    int i;
    double *p;

    p = coef;
    ans = *p++;
    i = N;
    
    do
      ans = ans * x  +  *p++;
    while( --i );
    
    return ans;
  }
  
  double DiLog(double x)
  {
    static double cof_A[8] = {
      4.65128586073990045278E-5,
      7.31589045238094711071E-3,
      1.33847639578309018650E-1,
      8.79691311754530315341E-1,
      2.71149851196553469920E0,
      4.25697156008121755724E0,
      3.29771340985225106936E0,
      1.00000000000000000126E0,
    };
    static double cof_B[8] = {
      6.90990488912553276999E-4,
      2.54043763932544379113E-2,
      2.82974860602568089943E-1,
      1.41172597751831069617E0,
      3.63800533345137075418E0,
      5.03278880143316990390E0,
      3.54771340985225096217E0,
      9.99999999999999998740E-1,
    };
    if( x >1. ) {
      return -DiLog(1./x)+M_PI*M_PI/3.-0.5*sqr(log(x));
    }
    x = 1.-x;
    double w, y, z;
    int flag;
    if( x == 1.0 )
      return( 0.0 );
    if( x == 0.0 )
      return( M_PI*M_PI/6.0 );
    
    flag = 0;
    
    if( x > 2.0 )
      {
	x = 1.0/x;
	flag |= 2;
      }
    
    if( x > 1.5 )
      {
	w = (1.0/x) - 1.0;
	flag |= 2;
      }
    
    else if( x < 0.5 )
      {
	w = -x;
	flag |= 1;
      }
    
    else
      w = x - 1.0;
    
    
    y = -w * polevl( w, cof_A, 7) / polevl( w, cof_B, 7 );
    
    if( flag & 1 )
	y = (M_PI * M_PI)/6.0  - log(x) * log(1.0-x) - y;
    
    if( flag & 2 )
      {
	z = log(x);
	y = -0.5 * z * z  -  y;
      }
    
    return y;

  }

}


