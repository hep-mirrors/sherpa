// ************************************************************************
// *  class Gauss_Integrator
// *   - provides routines for non MC Integration of a given 1D function by
// *      + Gauss-Legendre  W(x)=P(x) (General Polinominal Integrator)
// *      + Gauss-Chebyshev W(x)=(1-x^2)^-1/2 * P(x)
// *      + Gauss-Laguerre  W(x)=x^alf e^-x  * P(x)
// *      + Gauss-Hermite   W(x)=e^(-x^2) * P(x)
// *      + Gauss-Jacobi    W(x)=(1-x)^alf (1+x)^bet * P(x)
// ************************************************************************
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include <iostream>
#include <fstream>

using namespace ATOOLS;

int Gauss_Integrator::s_ngauleg=0;
int Gauss_Integrator::s_ngaulag=0;
int Gauss_Integrator::s_ngauherm=0;
int Gauss_Integrator::s_ngaujac=0;
Weight_Module* Gauss_Integrator::s_wlistroot=0;


Gauss_Integrator::Gauss_Integrator(Function_Base *func=0){
  m_numberabsc = 0;  // precisision (number of points) not jet choosen
  m_wlistact   = 0;    // so no precalculated weights or abscissas are available jet;
  m_func       = func;    // set function to be integrated
}

double Gauss_Integrator::Integrate(double x1, double x2, double prec, int mode, int nmax) {
  if (x1==x2) return 0.;

  int n = 64;
  if (n>nmax) n=nmax;
  double i2=0.,i1=1.;
  int err;
  for (;(n<=nmax)&&(dabs(1-i2/i1)>prec);n*=2) { 
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
      i1 = Jacobi(x1,x2,n,-0.5, -0.5);
      break;
    default: 
      i1 = Legendre(x1,x2,n);
    }
  }
  return i1;
}

// * Gauss_Legendre Integation ******************

void Gauss_Integrator::GauLeg(double * x, double * w, int n)
{
  const double EPS=3.0e-11;

  // calculate Gauss_Legendre weights and abscissas
  //  x[0] .. x[n-1]   abscissas
  //  w[0] .. w[n-1]   weights
  int m,j,i;
  double z1,z,pp,p3,p2,p1;
  
  m=(n+1)/2;
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
    x[i-1]=-z,x[n-i]=z;
    w[i-1]=2.0/((1.0-z*z)*pp*pp);
    w[n-i]=w[i-1];
  }
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
    if (n>s_ngauleg) {
      // generate new weights and abscissas
      m_wlistact    = new Weight_Module;
      m_wlistact->w = new double[n];
      m_wlistact->x = new double[n];
      m_wlistact->methode = 1;             // legendre
      m_wlistact->n       = n;             // number of points
      if (n>s_ngauleg) s_ngauleg=n;
      m_wlistact->next = s_wlistroot;      // put in list
      s_wlistroot=m_wlistact;
      GauLeg(m_wlistact->x, m_wlistact->w, n); 

    } 
    else {
      // use already calculated abscissas
      Weight_Module * listitem;
      listitem = s_wlistroot;
      m_wlistact = 0;
      for (;listitem!=0;listitem=listitem->next) 
	if ((listitem->n>=n) && (listitem->methode==1)) 
	  if ((m_wlistact==0) || (m_wlistact->n>listitem->n)) 
	    m_wlistact=listitem;
      if ((m_wlistact==0)||(m_wlistact->n>2*n)) { 
	// generate new weights and abscissas
	m_wlistact    = new Weight_Module;
	m_wlistact->w = new double[n];
	m_wlistact->x = new double[n];
	m_wlistact->methode = 1;             // legendre
	m_wlistact->n       = n;                 // number of points
	m_wlistact->next    = s_wlistroot;      // put in list
	s_wlistroot       = m_wlistact;
	if (n>s_ngauleg) s_ngauleg = n;
	GauLeg(m_wlistact->x, m_wlistact->w, n); 
      } ; 
    }
    // do the summation (the same for all gauss integrations);
    for (int i=0;i<n;i++) {
      double x=xm+xl*m_wlistact->x[i];
      sum+=m_wlistact->w[i]*((*m_func)( x ));
    }
    sum*=xl;
    return sum;
  }
}
 
double Gauss_Integrator::Chebyshev( double a, double b, double prec, int n_max, int &i_err ) 
{
// 	a;		lower boundary                            
// 	b;		upper boundary                            
// 	prec;		relative error                            
//      n_max;		maximum number of steps allowed           
//	*i_err;		= 0 for an alleged successful calculation 
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

  // initializing m, i_err, n, s0, c0, tsch and p. 
  intalt = intneu = 0.;

  h  = (b-a)/2.;
  m  = 0;
  n  = 1;
  s0 = 1.;
  c0 = 0.;
  t  = a + h;
  ch = (*m_func)( t );

  // computing the (2n+1) points quadrature formula. 
  // updating q, p, c1, s1, c0, s0, s and c. 
  while (m<5 || ( dabs( intneu-intalt ) > prec * dabs( intneu ) 
		   && m<n_max ) ) {
    intalt = intneu;
    c1 = c0;
    s1 = s0;
    c0 = sqrt( ( 1. + c1 ) * 0.5 );
    s0 = s1 / ( c0 + c0 );
    s  = s0;
    c  = c0;

    // computing f() at the new points. 
    for ( i = 1; i <= n; ){
      di = i;
      x = 1. + 0.21220659078919378103 * s * c * ( 3. + 2. * s * s )
	- di / ( n+1 );
      tm = a + h * ( -x + 1. );
      tp = a + h * ( x + 1. );
      ch = ch + ( (*m_func)( tm ) + (*m_func)( tp ) ) * pow( s, 4 );
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

  if (fabs(intneu - intalt)>prec * fabs(intneu))
    i_err = 1;
  else
    i_err = 0;
  
  ch = 16. * ch / ( 3. * (n+1) );
  ch = h * ch;
  
  return ch;
}


// ******* Gauss Jacobi Interation *************************************
void Gauss_Integrator::GauJac(double * x, double * w, int n, double alf, double bet)
{
  const double EPS   = 3.0e-14;
  const int    MAXIT = 10;
  // calculate Gauss Jacobi weights and abscissas
  int i,its,j;
  double alfbet,an,bn,r1,r2,r3;
  double a,b,c,p1,p2,p3,pp,temp,z=1.,z1;

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
    x[i-1]=z;
    w[i-1]=exp(Gammln(alf+n)+Gammln(bet+n)-Gammln(n+1.0)-
	     Gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2)
             *pow((1.-z),-alf)*pow((1.+z),-bet);
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software VsXz&v%120(9p+45$j3D. */

double Gauss_Integrator::Jacobi(double x1, double x2,int n=8, 
                               double alf= -0.5f, double bet= -0.5f)
{
  double sum=0.0;
  double xm=0.5*(x2+x1);
  double xl=0.5*(x2-x1);
  if (n>s_ngaujac) {
    // generate new weights and abscissas
    m_wlistact    = new Weight_Module;
    m_wlistact->w = new double[n];
    m_wlistact->x = new double[n];
    m_wlistact->methode = 5;             // Jacobi
    m_wlistact->n       = n;                 // number of points
    if (n>s_ngaujac) s_ngaujac = n;
    m_wlistact->next = s_wlistroot;      // put in list
    s_wlistroot      = m_wlistact;
    GauJac(m_wlistact->x, m_wlistact->w, n, alf, bet); 
  } 
  else {
    // use already calculated abscissas
    Weight_Module * listitem;
    listitem = s_wlistroot;
    m_wlistact = 0;
    for (;listitem!=0;listitem=listitem->next) 
      if ((listitem->n>=n)&&(listitem->methode==5)) 
	if ((m_wlistact==0)||(m_wlistact->n>listitem->n)) 
	  m_wlistact=listitem;
    if ((m_wlistact==0)||(m_wlistact->n>2*n)) { 
      // generate new weights and abscissas
      m_wlistact    = new Weight_Module;
      m_wlistact->w = new double[n];
      m_wlistact->x = new double[n];
      m_wlistact->methode = 5;             // jacobi
      m_wlistact->n       = n;             // number of points
      m_wlistact->next    = s_wlistroot;   // put in list
      s_wlistroot         = m_wlistact;
      if (n>s_ngaujac) s_ngaujac = n;
      GauJac(m_wlistact->x, m_wlistact->w, n, alf, bet); 
    }
  }
  
  // do the summation (the same for all gauss integrations);
  for (int i=0;i<n;i++) {
    double x=xm+xl*m_wlistact->x[i];
    sum+=m_wlistact->w[i]*((*m_func)( x ));
  }
  sum*=xl;
  return sum;
}









