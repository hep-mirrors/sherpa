#include "Traces.H"


// elementary matrices

double HADRONS::G( int mu, int nu )
{
  if( mu != nu ) return 0.;
  if( mu==0 ) return 1.;
  return -1.;
}

double HADRONS::E( short a, short b, short c, short d )
{
  if( a==b || a==c || a==d || b==c || b==d || c==d ) return 0.;
  short ind[4] = {a,b,c,d};
  short help;
  double ret=1.;
  // sort
  for( short r=1; r<4; r++ )
	for( short s=0; s<r; s++ )
	{
//	  if( ind[s]==ind[r] && r!=s ) return 0;
	  if( ind[s]>ind[r] ) {
		  help = ind[r];
		  ind[r] = ind[s];
		  ind[s] = help;
		  ret *= -1.;
	  }
	}
  return ret;
}
  
double HADRONS::Delta( int mu, int nu )
{
  if( mu==nu ) return 1.;
  return 0.;
}

// Traces with 2 gammas

double HADRONS::Trace( int mu, int nu )
{
  return 4.*G(mu,nu);
}

double HADRONS::Trace( Vec4D a, int B )
{
  return 4.*a[B];
}

double HADRONS::Trace( Vec4D a, Vec4D b )
{
  double dot = a*b;
  return 4.*dot;
}

// Traces with 4 gammas

double HADRONS::Trace( int mu, int nu, int rho, int sigma )
{
  return( 4.*(G(mu,nu)*G(rho,sigma)+G(nu,rho)*G(mu,sigma)-G(mu,rho)*G(nu,sigma)) );
}

double HADRONS::Trace( Vec4D a, int B, int C, int D )
{
  double term = 0.;
//  for( short A=0; A<4; A++ ) term += a[A]*Trace(A,B,C,D)*G(A,A);
  term = a[B]*G(C,D) + a[D]*G(B,C) - a[C]*G(B,D);
  return 4.*term;
}

double HADRONS::Trace( Vec4D a, Vec4D b, int C, int D )
{
  double term = 0.;
//  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) 
//	term += a[A]*b[B]*Trace(A,B,C,D)*G(A,A)*G(B,B);
  term = (a*b)*G(C,D) + a[D]*b[C] - a[C]*b[D];
  return 4.*term;
}

double HADRONS::Trace( Vec4D a, int C, Vec4D b, int D )
{
  double term = 0.;
//  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) 
//	term += a[A]*b[B]*Trace(A,C,B,D)*G(A,A)*G(B,B);
  term = a[C]*b[D] + a[D]*b[C] - (a*b)*G(C,D);
  return 4.*term;
}

double HADRONS::Trace( Vec4D a, Vec4D b, Vec4D c, int D )
{
  double term = 0.;
//  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) for( short C=0; C<4; C++ ) 
//	term += a[A]*b[B]*c[C]*Trace(A,B,C,D)*G(A,A)*G(B,B)*G(C,C);
  term = (a*b)*c[D] + (b*c)*a[D] - (a*c)*b[D];
  return 4.*term;
}

double HADRONS::Trace( Vec4D a, Vec4D b, Vec4D c, Vec4D d )
{
  double term = 0.;
  term = (a*b)*(c*d) + (b*c)*(a*d) - (a*c)*(b*d); 
  return 4.*term;
}

// Traces with 4 gammas and gamma5

Complex HADRONS::Trace5( int mu, int nu, int rho, int sigma )
{
  return Complex( 0., 4.*E(mu,nu,rho,sigma) );
}

Complex HADRONS::Trace5( Vec4D a, int B, int C, int D )
{
  if( B==C || B==D || C==D ) return Complex(0.,0.);
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) {
	if( A!=B && A!=C && A!=D ) term += a[A]*Trace5(A,B,C,D)*G(A,A);
  }
  return term;
}

Complex HADRONS::Trace5( Vec4D a, Vec4D b, int C, int D )
{
  if( C==D ) return Complex(0.,0.);
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) {
	if( A!=B && A!=C && A!=D && B!=C && B!=D ) term += a[A]*b[B]*Trace5(A,B,C,D)*G(A,A)*G(B,B);
  }
  return term;
}

Complex HADRONS::Trace5( Vec4D a, int C, Vec4D b, int D )
{
  if( C==D ) return Complex(0.,0.);
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) {
	if( A!=B && A!=C && A!=D && B!=C && B!=D ) term += a[A]*b[B]*Trace5(A,C,B,D)*G(A,A)*G(B,B);
  }
  return term;
}

Complex HADRONS::Trace5( Vec4D a, Vec4D b, Vec4D c, int D )
{
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) for( short C=0; C<4; C++ ) {
	if( A!=B && A!=C && A!=D && B!=C && B!=D && C!=D )
	  term += a[A]*b[B]*c[C]*Trace5(A,B,C,D)*G(A,A)*G(B,B)*G(C,C);
  }
  return term;
}

Complex HADRONS::Trace5( Vec4D a, Vec4D b, Vec4D c, Vec4D d )
{
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) for( short C=0; C<4; C++ ) for( short D=0; D<4; D++ ) {
	if( A!=B && A!=C && A!=D && B!=C && B!=D && C!=D )
	  term += a[A]*b[B]*c[C]*d[D]*Trace5(A,B,C,D)*G(A,A)*G(B,B)*G(C,C)*G(D,D);
  }
  return term;
}


// Traces with 6 gammas

double HADRONS::Trace( int mu, int nu, int rho, int sigma, int alpha, int beta )
{
  double ret=0.;
  double gf = G(mu,mu);
  if( mu==nu ) 		ret += Trace(rho,sigma,alpha,beta); 
  if( mu==rho )		ret -= Trace(nu,sigma,alpha,beta); 
  if( mu==sigma )	ret += Trace(nu,rho,alpha,beta); 
  if( mu==alpha )	ret -= Trace(nu,rho,sigma,beta); 
  if( mu==beta )	ret += Trace(nu,rho,sigma,alpha); 
  return ret*gf;
}

double HADRONS::Trace( Vec4D a, int mu, int nu, Vec4D b, int rho, int sigma )
{
  double term(0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) {
	term += a[A]*b[B]*Trace(A,mu,nu,B,rho,sigma)*G(A,A)*G(B,B);
  }
  return term;
}

// Traces with 6 gammas and gamma5

Complex HADRONS::Trace5( int mu, int nu, int rho, int sigma, int alpha, int beta ) 
{
  Complex realpart(0.,0.);
  if( alpha == beta && mu!=nu && mu!=rho && mu!=sigma && nu!=rho && nu!=sigma && rho!=sigma ) 
	realpart = Trace5(mu,nu,rho,sigma);
  double imagpart(0.);
  for( int A=0; A<4; A++ ) for( int B=0; B<4; B++ ) {
	if( A!=B && A!=alpha && A!=beta && B!=alpha && B!=beta && alpha!=beta )
	  imagpart += G(A,A)*G(B,B)*E(A,B,alpha,beta)*Trace(mu,nu,rho,sigma,A,B);
  }
  return( realpart + Complex(0.,imagpart/2.) );
}

Complex HADRONS::Trace5( Vec4D a, int mu, int nu, Vec4D b, int rho, int sigma )
{
  Complex term(0.,0.);
  for( short A=0; A<4; A++ ) for( short B=0; B<4; B++ ) {
	term += a[A]*b[B]*Trace5(A,mu,nu,B,rho,sigma)*G(A,A)*G(B,B);
  }
  return term;
}


// Propagator structure function

double HADRONS::Zeta( Vec4D a, Vec4D b, Vec4D c, double MP2 )
{
  double d1 = a*b;
  double d2 = a*c;
  double d3 = b*c;
  return d1 - d2*d3/MP2;
}

double HADRONS::Zeta( int mu, int nu, Vec4D c, double MP2 )
{
  return G(mu,nu)-c[mu]*c[nu]/MP2;
}
