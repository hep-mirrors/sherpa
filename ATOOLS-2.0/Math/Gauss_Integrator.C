// ************************************************************************
// *  class Gauss_Integrator
// *   - provides routines for non MC Integration of a given 1D function by
// *      + Gauss-Legendre  W(x)=P(x) (General Polinominal Integrator)
// *      + Gauss-Chebyshev W(x)=(1-x^2)^-1/2 * P(x)
// *      + Gauss-Laguerre  W(x)=x^alf e^-x  * P(x)
// *      + Gauss-Hermite   W(x)=e^(-x^2) * P(x)
// *      + Gauss-Jacobi    W(x)=(1-x)^alf (1+x)^bet * P(x)
// ************************************************************************
#include "Gauss_Integrator.H"
#include "MathTools.H"
#include <iostream>
#include <fstream>

using namespace AMATOOLS;
using std::cout;
using std::cerr;
using std::endl;

int Gauss_Integrator::n_gauleg=0;
int Gauss_Integrator::n_gaulag=0;
int Gauss_Integrator::n_gauherm=0;
int Gauss_Integrator::n_gaujac=0;
weightlist* Gauss_Integrator::wlistroot=0;


Gauss_Integrator::Gauss_Integrator(Function_Base *_func=0){
  mode=0;        // Methode not jet choosen
  numberabsc=0;  // precisision (number of points) not jet choosen
  wlistact=0;    // so no precalculated weights or abscissas are available jet;
  func=_func;    // set function to be integrated
}

double Gauss_Integrator::Integrate(double x1, double x2, double prec, int mode, int nmax) {
  if (x1==x2) return 0.;

  int n = 64;
  //  double guess=-0.25*log(prec)/log(2);
  //  n=n* int(pow(2.,int(guess))+0.5);
  //  cout<<" guess n="<<n<<endl;
  if (n>nmax) {
    n=nmax;
    cout<<" reduced to n="<<n<<endl;
  }
  double i2=0.,i1=1.;
  int err;
  for (n;(n<=nmax)&&(dabs(1-i2/i1)>prec);n*=2) { 
    i2=i1;
    switch (mode) {
    case 1 :
      i1=Legendre(x1,x2,n);
      break;
    case 2 : 
      return Chebyshev(x1,x2,prec,n*4,err);
      // case 3 : return Laguerre(n);
      // case 4 : return Hermite(n);
    case 5 : 
      i1=Jacobi(x1,x2,n,-0.5, -0.5);
      break;
    default: 
      i1= Legendre(x1,x2,n);
    }
    //    cout<<"Guess the "<<n<<"th="<<i1<<" ("<<(dabs(i1-i2)/i1)<<")"<<endl;
    //cout<<"LegG("<<n<<".)="<<myintegrator.Legendre(s1,s2,n)<<endl;
  }
  return i1;
}

// * Gauss_Legendre Integation ******************

void Gauss_Integrator::gauleg(double x[], double w[], int n)
{
  const double EPS=3.0e-11;

  // calculate Gauss_Legendre weights and abscissas
  //  x[0] .. x[n-1]   abscissas
  //  w[0] .. w[n-1]   weights
  int m,j,i;
  double z1,z,pp,p3,p2,p1;
  
  m=(n+1)/2;
  //	xm=0.5*(x2+x1);
  //      xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
    z=cos(3.141592654*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    // x[i]=xm-xl*z;
    // x[n+1-i]=xm+xl*z;
    x[i-1]=-z,x[n-i]=z;
    // w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[i-1]=2.0/((1.0-z*z)*pp*pp);
    w[n-i]=w[i-1];
  }
  /*
  if ((n==3073)||(n==13)) {
    cout<<"writing legendre weights and abscissas to file.."<<endl;
    std::ofstream to;
    if (n==3072) to.open("gauleg.3072.dat");
    else to.open("gauleg.12.dat");
    to<<"fortran output:"<<endl<<endl;
    to<<"      DATA wg"<<endl;
    to<<"     $/";
    to.precision(15);
    //    to.setf(std::ios_base::scientific,std::ios_base::floatfield);
    to.setf(std::ios::scientific,std::ios::floatfield);
    for (int j=0; j<n-1 ;++j) {
      if (j%3 ==2) to<<w[j]<<","<<endl<<"     $ ";
      else to<<w[j]<<", ";
    }
    to<<w[n-1]<<"/"<<endl;

    to<<"      DATA xx"<<endl;
    to<<"     $/";
    for (int j=0; j<n-1 ;++j) {
      if (j%3 ==2) to<<x[j]<<","<<endl<<"     $ ";
      else to<<x[j]<<", ";
    }
    to<<x[n-1]<<"/"<<endl;
  }
  */
}
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

double Gauss_Integrator::Legendre(double x1, double x2,int n=8)
{
  if (n>32) {
    double xm=0.5*(x1 + x2);
    int nn = n/2;
    return Legendre(x1,xm,nn)+Legendre(xm,x2,nn);
  }
  else {

    double sum=0.0;
    double xm=0.5*(x2+x1);
    double xl=0.5*(x2-x1);
    // cout<<"already "<<n_gauleg<<" points calculated\n";
    if (n>n_gauleg) {
      // generate new weights and abscissas
      wlistact = new weightlist;
      wlistact->w =new double[n];
      wlistact->x =new double[n];
      wlistact->methode=1;             // legendre
      wlistact->n = n;                 // number of points
      if (n>n_gauleg) n_gauleg=n;
      wlistact->next = wlistroot;      // put in list
      wlistroot=wlistact;
      cout<<"calculate "<<n<<" legendre weights and abscissas..."<<endl;
      gauleg(wlistact->x, wlistact->w, n); 

    } else {
      // use already calculated abscissas
      // cout<<"looking for weights and absiscissas ...\n";
      weightlist * listitem;
      listitem = wlistroot;
      wlistact = 0;
      for (;listitem!=0;listitem=listitem->next) 
	if ((listitem->n>=n)&&(listitem->methode==1)) 
	  if ((wlistact==0)||(wlistact->n>listitem->n)) 
	    wlistact=listitem;
      if ((wlistact==0)||(wlistact->n>2*n)) { 
	cout<<" * no aproriate weights found!\n";
	// generate new weights and abscissas
	wlistact = new weightlist;
	wlistact->w =new double[n];
	wlistact->x =new double[n];
	wlistact->methode=1;             // legendre
	wlistact->n = n;                 // number of points
	wlistact->next = wlistroot;      // put in list
	wlistroot=wlistact;
	if (n>n_gauleg) n_gauleg=n;
	cout<<"calculate "<<n<<" new legendre weights and abscissas...\n";
	gauleg(wlistact->x, wlistact->w, n); 
      } ; //else cout<<"...found!\n";
    }
    // do the summation (the same for all gauss integrations);
    if (n!=(wlistact->n)) cout<<" increase n to "<<(n=(wlistact->n))<<" from "<<n<<endl;
    for (int i=0;i<n;i++) {
      double x=xm+xl*wlistact->x[i];
      sum+=wlistact->w[i]*((*func)( x ));
      //      cout<<" "<<i<<": x="<<x<<" w="<<wlistact->w[i]
      //        	<<" f="<<(*func)(x)<<" sum="<<sum<<endl;
    }
    sum*=xl;
    return sum;
  }
}

/*  great integration routine, from Ch. Hofmann.
 *  Changed to integrate classes derived from dfunc.       
 *                                      M. W. Beinker, 27.02.1996
 */

//#include "df.H" 
//extern "C" { 
//#include <math.h> 
//} 
 
double Gauss_Integrator::Chebyshev( double a, double b, double prec, int N_Max, int &I_Err ) 
{
// 	a;		lower boundary                            
// 	b;		upper boundary                            
// 	prec;		relative error                            
//      N_Max;		maximum number of steps allowed           
//	*I_Err;		= 0 for an alleged successful calculation 
//			= 1 otherwise                             
//      f;              dfunc (integrated function)               
  double 	ch;		// value of integral 

//      ref. 1: j. m. perez-jorda, e. san-fabian, f. moscardo,    
//              comp. phys. comm. 70 (1992) 271	                  
  int m, n, i; 
  double di;
  double s, s0, s1, c, c0, c1, tm, tp, x, t, h;
  double intalt;
  double intneu;


  // initializing m, I_Err, n, s0, c0, tsch and p. 

  intalt=intneu = 0.;

  h = (b-a)/2.;
  m = 0;
  n = 1;
  s0 = 1.;
  c0 = 0.;
  t = a + h;
  ch = (*func)( t );
  //p = ch;

  // computing the (2n+1) points quadrature formula. 
  // updating q, p, c1, s1, c0, s0, s and c. 

  while( m < 5 || ( dabs( intneu-intalt ) > prec * dabs( intneu ) 
		   && m < N_Max ) ) {

    intalt = intneu;
    //    q = p + p;
    //p = ch + ch;
    c1 = c0;
    s1 = s0;
    c0 = sqrt( ( 1. + c1 ) * 0.5 );
    s0 = s1 / ( c0 + c0 );
    s = s0;
    c = c0;

    // computing f() at the new points. 

    for ( i = 1; i <= n; ){
      di = i;
      x = 1. + 0.21220659078919378103 * s * c * ( 3. + 2. * s * s )
	- di / ( n+1 );
      tm = a + h * ( -x + 1. );
      tp = a + h * ( x + 1. );
      ch = ch + ( (*func)( tm ) + (*func)( tp ) ) * pow( s, 4 );
      x = s;
      s = s * c1 + c * s1;
      c = c * c1 - x * s1;
      i = i + 2;
    }

    // replacing n by 2n+1.

    m = m + 1;
    n = n + n + 1;
    
    intneu = ch / ( n+1 );
  }

  // test for successfullness and integral final value. 
  cout<<" used "<<n<<" points in Chebyshev"<<endl;

  if ( fabs( intneu - intalt ) > prec * fabs( intneu ) )
    I_Err = 1;
  else
    I_Err = 0;
  
  ch = 16. * ch / ( 3. * (n+1) );
  ch = h * ch;
  
  return ch;
}


// ******* Gauss Jacobi Interation *************************************
void Gauss_Integrator::gaujac(double x[], double w[], int n, double alf, double bet)
{
  const double EPS   = 3.0e-14;
  const int    MAXIT = 10;
  // calculate Gauss Jacobi weights and abscissas
  //    double gammln(double xx);
  //	void nrerror(char error_text[]);
  int i,its,j;
  double alfbet,an,bn,r1,r2,r3;
  double a,b,c,p1,p2,p3,pp,temp,z,z1;

  for (i=1;i<=n;i++) {
    if (i == 1) {
      an=alf/n;
      bn=bet/n;
      r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
      r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
      z=1.0-r1/r2;
    } else if (i == 2) {
      r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
      r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
      r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
      z -= (1.0-z)*r1*r2*r3;
    } else if (i == 3) {
      r1=(1.67+0.28*alf)/(1.0+0.37*alf);
      r2=1.0+0.22*(n-8.0)/n;
      r3=1.0+8.0*bet/((6.28+bet)*n*n);
      z -= (x[0]-z)*r1*r2*r3;
    } else if (i == n-1) {
      r1=(1.0+0.235*bet)/(0.766+0.119*bet);
      r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
      r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
      z += (z-x[n-4])*r1*r2*r3;
    } else if (i == n) {
      r1=(1.0+0.37*bet)/(1.67+0.28*bet);
      r2=1.0/(1.0+0.22*(n-8.0)/n);
      r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
      z += (z-x[n-3])*r1*r2*r3;
    } else {
      z=3.0*x[i-2]-3.0*x[i-3]+x[i-4];
    }
    alfbet=alf+bet;
    for (its=1;its<=MAXIT;its++) {
      temp=2.0+alfbet;
      p1=(alf-bet+temp*z)/2.0;
      p2=1.0;
      for (j=2;j<=n;j++) {
	p3=p2;
	p2=p1;
	temp=2*j+alfbet;
	a=2*j*(j+alfbet)*(temp-2.0);
	b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
	c=2.0*(j-1+alf)*(j-1+bet)*temp;
	p1=(b*p2-c*p3)/a;
      }
      pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS) break;
    }
    if (its > MAXIT) cout<<"too many iterations in Gauss_Integrator::gaujac";
    x[i-1]=z;
    w[i-1]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
	     gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2)
             *pow((1.-z),-alf)*pow((1.+z),-bet);
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */


double Gauss_Integrator::gammln(double xx)
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

double Gauss_Integrator::Jacobi(double x1, double x2,int n=8, 
                               double alf= -0.5f, double bet= -0.5f)
{
  double sum=0.0;
  double xm=0.5*(x2+x1);
  double xl=0.5*(x2-x1);
  if (n>n_gaujac) {
    // generate new weights and abscissas
    wlistact = new weightlist;
    wlistact->w =new double[n];
    wlistact->x =new double[n];
    wlistact->methode=5;             // Jacobi
    wlistact->n = n;                 // number of points
    if (n>n_gaujac) n_gaujac=n;
    wlistact->next = wlistroot;      // put in list
    wlistroot=wlistact;
    cout<<"calculate "<<n<<" jacobi weights and abscissas..."<<endl;
    gaujac(wlistact->x, wlistact->w, n, alf, bet); 
  } else {
    // use already calculated abscissas
    // cout<<"looking for weights and absiscissas ...\n";
    weightlist * listitem;
    listitem = wlistroot;
    wlistact = 0;
    for (;listitem!=0;listitem=listitem->next) 
      if ((listitem->n>=n)&&(listitem->methode==5)) 
	if ((wlistact==0)||(wlistact->n>listitem->n)) 
	  wlistact=listitem;
    if ((wlistact==0)||(wlistact->n>2*n)) { 
      cout<<" * no aproriate weights found!\n";
      // generate new weights and abscissas
      wlistact = new weightlist;
      wlistact->w =new double[n];
      wlistact->x =new double[n];
      wlistact->methode=5;             // jacobi
      wlistact->n = n;                 // number of points
      wlistact->next = wlistroot;      // put in list
      wlistroot=wlistact;
      if (n>n_gaujac) n_gaujac=n;
      cout<<"calculate "<<n<<" new legendre weights and abscissas..."<<endl;
      gaujac(wlistact->x, wlistact->w, n, alf, bet); 
    }  //else cout<<"...found!\n";
  }
  
  // do the summation (the same for all gauss integrations);
  for (int i=0;i<n;i++) {
    double x=xm+xl*wlistact->x[i];
    sum+=wlistact->w[i]*((*func)( x ));

    //  cout<<" "<<i<<": x="<<x<<" w="<<wlistact->w[i]
    //	<<" f="<<(*func)(x)<<" sum="<<sum<<endl;
  }
  sum*=xl;
  return sum;
}









