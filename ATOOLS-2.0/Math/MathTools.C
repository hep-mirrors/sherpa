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

}


