#include "ATOOLS/Math/MathTools.H"

#include "ATOOLS/Org/Message.H"

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

  Complex DiLog(const Complex& x) {
  // only calculates for real arguments
  // -> be careful with discontinuity along pos. real axis
    if (imag(x) == 0.) {
      std::cout<<"use real dilog ..."<<std::endl;
      return DiLog(real(x));
    }
    else {
      double precission = 1E-20;
      // needs to be larger than 0.5
      double radius     = 0.8;
      //! ref.: Mathematica(TM)
      //! and
      //! ref.: D.Maitre; "Extension of HPL to complex arguments";
      //!       hep-ph/0703052
      //! and
      //! ref.: C.Osacar, J.Palacian, M.Palacios; "Numerical Evaluation
      //!       of the Dilogarithm of Complex Argument"; Celestial Mechanics
      //!       and Dynamical Astronomy 62: 93-98, 1995
      //! maybe alternatively implement continued-fraction expansion of
      //! ref.: D.Cvijovic, J.Klinowski; "Continued-Fraction Expansions for
      //!       the Riemann Zeta Function and Polylogarithms"; Proceedings
      //!       of the American Mathematical Society, Volume 125, Number 9,
      //!       September 1997, Pages 2543-2550
      if (abs(x) < radius) {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region I ..."<<std::endl;
        // DiLog(x) = sum_{1}^{infty} x^i/i^2
        // 50 terms give rel. error of < 6.7E-19 using 19 dec. places
        Complex prod(x);
        Complex next(prod);
        Complex sum(next);
        unsigned int i=1;
        // to do: smart up precission test to include less terms when possible
        while /*((real(next)/real(sum) > precission) ||
               (imag(next)/imag(sum) > precission))*/(i < 100) {
          ++i;
          prod *= x;
          next  = prod/((double)(i*i));
          sum  += next;
        }
        std::cout<<"number of iterations: "<<i<<std::endl;
        return sum;
        // alternatively DiLog(x) = 4x2/(1+4x+x2)*sum_{1}^{infty} x^i/(i(i+1)(i+2))^2
        //                          + 4x + 23/4 x2 + 3(1-x2)ln(1-x)
        // exclude vicinity of -2+sqrt(3)=-0.2....
//         Complex prod(x);
//         Complex next(prod/24.);
//         Complex sum(next);
//         Complex factor(4.*x*x/(1.+4.*x+x*x));
//         double convfactor(1E-12/abs(factor));
//         unsigned int i=1;
//         while (abs(next/sum) > convfactor) {
//           ++i;
//           prod *= x;
//           next  = prod/sqr((double)(i*(i+1)*(i+2)));
//           sum  += next;
//         }
//         return factor*sum + 4.*x + 23./4.*x*x + 3.*(1.-x*x)*log(1.-x);
      }
      else if (abs(x) > 1./radius) {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region II ..."<<std::endl;
        // DiLog(x) = -DiLog(1/x) - 1/2*ln2(-x) - pi2/6
        return -DiLog(1./x) - 0.5*sqr(log(-x)) - sqr(M_PI)/6.;
      }
      else if (abs(1.-x) < radius) {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region III ..."<<std::endl;
        // DiLog(x) = -DiLog(1-x) - ln x ln(1-x) + pi2/6
        return -DiLog(1.-x) - log(x)*log(1.-x) + sqr(M_PI)/6.;
      }
      else if (abs(x+1.) < radius) {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region IV ..."<<std::endl;
        // DiLog(x) = -DiLog(-x/(1-x)) - 1/2 ln2(1-x)
//         return -DiLog(-x/(1.-x)) - 0.5*sqr(log(1.-x));
        // DiLog(x) = -DiLog(1-x)-log(x)*log(1-x)+pi2/6;
        return -DiLog(1.-x)-log(x)*log(1.-x) + sqr(M_PI)/6.;
      }
      else if (abs(x-Complex(0.,1.)) < 0.7*radius) {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region V ..."<<std::endl;
        // expansion around (0,i), coeffients taken from Mathematica(TM)
        Complex A[35] = {
          Complex(-0.20561675835602830456,      0.91596559417721901505),      // y^0
          Complex( 0.78539816339744830962,      0.34657359027997265471),      // y^1
          Complex( 0.07671320486001367265,      0.14269908169872415481),      // y^2
          Complex(-0.011799387799149436539,     0.051142136573342448430 ),    // y^3
          Complex(-0.017523269096673502990,     0.011983792483971255929),     // y^4

          Complex(-0.0095870339871770047435,   -0.0015186152773388023916 ),   // y^5
          Complex(-0.0029011539355509980070,   -0.0038225283226475039529),    // y^6
          Complex( 0.0003002623717930986263,   -0.0024867033733294268631),    // y^7
          Complex( 0.00105979402309181993382,  -0.00085334185325246727339),   // y^8
          Complex( 0.00075852609177997090968,   0.00007398357608161771895),   // y^9

          Complex( 0.00028063700374876627517,   0.00033545126037975159649),   // y^10
          Complex(-0.00002086478216341054226,   0.00025512454886251479561),   // y^11
          Complex(-0.00011549295766942644143,   0.00009924516180475245747),   // y^12
          Complex(-0.000091610918589002268437, -6.448627592291074140E-6),     // y^13
          Complex(-0.000036937812840125321265, -0.000042141457371106502010),  // y^14

          Complex( 2.129645927318449495E-6,    -0.000034475291984116966514),  // y^15
          Complex( 0.000016044544568442989441, -0.000014279498609805620265),  // y^16
          Complex( 0.0000134395281033464661318, 7.395345938286959441E-7 ),    // y^17
          Complex( 5.6843120992925061182E-6,    6.3101262152520546147E-6),    // y^18
          Complex(-2.671224963206599274E-7,     5.3851377782771110593E-6 ),   // y^19

          Complex(-2.5459795735737818222E-6,    2.3161349442848467532E-6),    // y^20
          Complex(-2.2058428040808064316E-6,   -9.95936414988398306E-8),      // y^21
          Complex(-9.618191647164753998E-7,    -1.0486913995663109012E-6),    // y^22
          Complex( 3.811341381441596081E-8,    -9.2000094016358516500E-7),    // y^23
          Complex( 4.3938382671474013095E-7,   -4.0575838603654702307E-7),    // y^24

          Complex( 3.8952805059508514214E-7,    1.490743197948385905E-8),     // y^25
          Complex( 1.7346641155818859707E-7,    1.8674572172604340591E-7),    // y^26
          Complex(-5.93987946838367862E-9,      1.6704172964862605644E-7),    // y^27
          Complex(-8.0341619925725390401E-8,    7.5006592748079759635E-8),    // y^28
          Complex(-7.2420158515387354131E-8,   -2.404770654887569057E-9),     // y^29

          Complex(-3.2753064372689143015E-8,   -3.4928477225793982556E-8),    // y^30
          Complex( 9.87152019585574517E-10,    -3.1696513909054009369E-8),    // y^31
          Complex( 1.5324154036392039318E-8,   -1.4425540294030506945E-8),    // y^32
          Complex( 1.39884027093629158254E-8,   4.101748776187957026E-10),    // y^33
          Complex( 6.4017059621418284134E-9,    6.7771622274921705293E-9)     // y^34
        };
        Complex y(x-Complex(0.,1.));
        Complex prod(y);
        Complex next(prod*A[1]);
        Complex sum(A[0]);
        unsigned int i=1;
        while (((real(next)/real(sum) > precission) ||
                (imag(next)/imag(sum) > precission)) && (i<35)) {
          ++i;
          sum  += next;
          prod *= y;
          next  = prod*A[i];
        }
        std::cout<<"number of iterations: "<<i<<std::endl;
        return sum;
      }
      else {
        std::cout<<abs(x)<<std::endl;
        std::cout<<"region VI ..."<<std::endl;
        // DiLog(x) = -DiLog(1/x) - 1/2*ln2(-x) - pi2/6
        return -DiLog(1./x) - 0.5*sqr(log(-x)) - sqr(M_PI)/6.;
        // DiLog(x) = -DiLog(-x) + 1/2*DiLog(x2)
	// return -DiLog(-x) + 0.5*DiLog(x*x);
      }
    }
    std::cout<<"no suitable region found ..."<<std::endl;
    return 0.;
  }

  int Factorial(const int n)
  {
    if (n < 0) return 0;
    if (n < 2) return 1;
    int res = 2;
    for (int i=3;i<=n;i++) res *= i;
    return res;
  }

  double ExpIntegral(int n, double x) {
    // Implementation close to chapter 6.3 in
    // "NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
    // (ISBN 0-521-43108-5)

    // configuration
    const int MaximumIterations = 100;
    const double EulerMascheroni = 0.577215664901532860606512;

    const double eps = 1.e-4; // relative error

    // Check arguments
    if (n < 0 || x < 0.0 || ((x==0.0) && ((n==0) || (n==1)))) {
      msg_Error() << "Bad arguments in E_n(x)" << std::endl;
    }
    else {
      // Handle special cases
      if (n == 0) return exp(-x)/x; // n=0
      else if (fabs(x) < 1.e-10) { // x=0
	return 1.0/(n-1);
      }
      // Normal Calculation
      else if (x > 1.) { // continued fraction
	double b=x+n;
	double c=1.e30; // 1./DBL_MIN (<cfloats>)
	double d=1./b;
	double h=d;
	for (int i=1; i <= MaximumIterations; i++) {
	  double a = -i*(n-1+i);
	  b+=2.;
	  d=1./(a*d + b);
	  c=b + a/c;
	  double del=c*d;
	  h *= del;

	  if (fabs(del-1.0) < eps) {
	    return h*exp(-x);
	  }
	}
	msg_Error() << "Continued fraction failed in ExpIntegral()! x=" << x
		    << std::endl;
      }
      else { // series evaluation
	double result = (n-1) ? 1./(n-1) : -log(x)-EulerMascheroni;

	double fact=1.;
	for (int i=1; i <= MaximumIterations; i++) {
	  fact *= -x/i;
	  double del = 0.;
	  if (i != (n - 1)) del = -fact/(i-n+1);
	  else {
	    double psi = -EulerMascheroni;
	    for (int j=1; j < (n - 1); j++) {
	      psi += 1./j;
	    }
	    del=fact*(-log(x)+psi);
	  }
	  result += del;
	  if (fabs(del) < fabs(result)*eps) {
	    return result;
	  }
	}
	msg_Error() << "Series failed in ExpIntegral()! x=" << x << std::endl;
      }
    }
    msg_Error() << "Exponential Integral Calculation failed! x=" << x
		<< std::endl;
    return 0.0;
  }

  double evaluate_polynomial (size_t size, double P[], double x)
{
	double value = 0;
	int deg = size-1;
	for (int i=deg; i>=0; i--) {
		value += P[i]*pow(x,i);
	}
	return value;
}

long double cyl_bessel_0 (long double x) {

  if (x>35) { //accurate to 1/1000
    long double ret = expl(-x)*sqrt(M_PI/(2*x)) * (1 - 1/(8*x)); // limiting expression
    //std::cout <<"Using limiting expression, bessel_0 = " << ret << std::endl; //debugging
    return ret;
  }

	static double P1[] = {
	2.4708152720399552679e+03,
	5.9169059852270512312e+03,
	4.6850901201934832188e+02,
	1.1999463724910714109e+01,
	1.3166052564989571850e-01,
	5.8599221412826100000e-04
	};

	static double Q1[] = {
	2.1312714303849120380e+04,
	-2.4994418972832303646e+02,
	1.0
	};

	static double P2[] = {
	-1.6128136304458193998e+06,
	-3.7333769444840079748e+05,
	-1.7984434409411765813e+04,
	-2.9501657892958843865e+02,
	-1.6414452837299064100e+00
	};

	static double Q2[] = {
	-1.6128136304458193998e+06,
	2.9865713163054025489e+04,
	-2.5064972445877992730e+02,
	1.0
	};

	static double P3[] = {
	1.1600249425076035558e+02,
	2.3444738764199315021e+03,
	1.8321525870183537725e+04,
	7.1557062783764037541e+04,
	1.5097646353289914539e+05,
	1.7398867902565686251e+05,
	1.0577068948034021957e+05,
	3.1075408980684392399e+04,
	3.6832589957340267940e+03,
	1.1394980557384778174e+02
	};

	static double Q3[] = {
	9.2556599177304839811e+01,
	1.8821890840982713696e+03,
	1.4847228371802360957e+04,
	5.8824616785857027752e+04,
	1.2689839587977598727e+05,
	1.5144644673520157801e+05,
	9.7418829762268075784e+04,
	3.1474655750295278825e+04,
	4.4329628889746408858e+03,
	2.0013443064949242491e+02,
	1.0
	};

	double factor, r, r1, r2;
  long double value;

	if (x < 0) return -1; //error
	if (x == 0) return -1; //error

	if (x <= 1) {
		double y = x*x;
		r1 = evaluate_polynomial(sizeof(P1)/sizeof(P1[0]),P1,y) / evaluate_polynomial(sizeof(Q1)/sizeof(Q1[0]),Q1,y);
		r2 = evaluate_polynomial(sizeof(P2)/sizeof(P2[0]),P2,y) / evaluate_polynomial(sizeof(Q2)/sizeof(Q2[0]),Q2,y);
		factor = log(x);
		value = r1 - factor * r2;
	}
	else {
		double y = 1/x;
		r = evaluate_polynomial(sizeof(P3)/sizeof(P3[0]),P3,y) / evaluate_polynomial(sizeof(Q3)/sizeof(Q3[0]),Q3,y);
		factor = exp(-x)/sqrt(x);
		value = factor * r;
	}
  //std::cout << "Using rational approximation, bessel_0 = " << value << std::endl; //debugging
	return value;
}

long double cyl_bessel_1 (long double x) {

  if (x>35) { //accurate to 1/1000
    long double ret = expl(-x)*sqrt(M_PI/(2*x)) * (1 + 3/(8*x)); // limiting expression
    //std::cout << "Using limiting expression, bessel_1 = " << ret << std::endl; //debugging
    return ret;
  }

	static double P1[] = {
	-2.2149374878243304548e+06,
	7.1938920065420586101e+05,
	1.7733324035147015630e+05,
	7.1885382604084798576e+03,
	9.9991373567429309922e+01,
	4.8127070456878442310e-01
	};

	static double Q1[] = {
	-2.2149374878243304548e+06,
	3.7264298672067697862e+04,
	-2.8143915754538725829e+02,
	1.0
	};

	static double P2[] = {
	0.0,
	-1.3531161492785421328e+06,
	-1.4758069205414222471e+05,
	-4.5051623763436087023e+03,
	-5.3103913335180275253e+01,
	-2.2795590826955002390e-01
	};

	static double Q2[] {
	-2.7062322985570842656e+06,
	4.3117653211351080007e+04,
	-3.0507151578787595807e+02,
	1.0
	};

	static double P3[] = {
	2.2196792496874548962e+00,
	4.4137176114230414036e+01,
	3.4122953486801312910e+02,
	1.3319486433183221990e+03,
	2.8590657697910288226e+03,
	3.4540675585544584407e+03,
	2.3123742209168871550e+03,
	8.1094256146537402173e+02,
	1.3182609918569941308e+02,
	7.5584584631176030810e+00,
	6.4257745859173138767e-02
	};

	static double Q3[] = {
	1.7710478032601086579e+00,
	3.4552228452758912848e+01,
	2.5951223655579051357e+02,
	9.6929165726802648634e+02,
	1.9448440788918006154e+03,
	2.1181000487171943810e+03,
	1.2082692316002348638e+03,
	3.3031020088765390854e+02,
	3.6001069306861518855e+01,
	1.0
	};

	double factor, r, r1, r2, y;
  long double value;

	if (x < 0) return -1; //error
	if (x == 0) return -1; //error

	if (x <= 1) {
		y = x*x;
		r1 = evaluate_polynomial(sizeof(P1)/sizeof(P1[0]),P1,y) / evaluate_polynomial(sizeof(Q1)/sizeof(Q1[0]),Q1,y);
		r2 = evaluate_polynomial(sizeof(P2)/sizeof(P2[0]),P2,y) / evaluate_polynomial(sizeof(Q2)/sizeof(Q2[0]),Q2,y);
		factor = log(x);
		value = (r1 + factor*r2)/x;
	}
	else {
		y = 1/x;
		r = evaluate_polynomial(sizeof(P3)/sizeof(P3[0]),P3,y) / evaluate_polynomial(sizeof(Q3)/sizeof(Q3[0]),Q3,y);
		factor = exp(-x)/sqrt(x);
		value = factor*r;
	}

  //std::cout << "Using rational approximation, bessel_1 = " << value << std::endl; //debugging
	return value;
}

long double cyl_bessel_2 (long double x) {
	return 2*1/x * cyl_bessel_1(x) + cyl_bessel_0(x);
}

}
