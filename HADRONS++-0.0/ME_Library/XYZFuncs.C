#include "Vector.H"
#include "Flavour.H"
#include "XYZFuncs.H"

using namespace HADRONS;

Complex XYZFunc::S( const int s, const int i, const int j )
{
  Complex A = eta[j]/eta[i];
  Complex ret(0.,0.);
  Vec4D pi, pj;
  pi = p_mom[i];
  pj = p_mom[j];
  ret  = Complex( 1.*(double)s*pi[1], SQRT_05*(pi[2]-pi[3]) )*A;
  ret -= Complex( 1.*(double)s*pj[1], SQRT_05*(pj[2]-pj[3]) )/A;
  return ret;
}  

Complex XYZFunc::Z( const int t1, const int t2, const int t3, const int t4,
			   const int hel_comb,
			   const Complex cR1, const Complex cL1,
			   const Complex cR2, const Complex cL2 )
{
  Complex z(0.,0.);
  switch( hel_comb ) {
  case 0	: z  = S(1,t3,t1)*S(-1,t4,t2)*cR1*cR2;
			  z -= mu[t1]*mu[t2]*eta[t3]*eta[t4]*cL1*cR2;
			  z -= mu[t3]*mu[t4]*eta[t1]*eta[t2]*cR1*cL2;
			  z *= -2.;
			  break;
  case 1	: z  = S(1,t4,t1)*mu[t3]*cL2;
			  z -= S(+1,t3,t1)*mu[t4]*cR2;
			  z *= -2.*eta[t2]*cR1;
			  break;
  case 2	: z  = S(-1,t2,t3)*mu[t4]*cL2;
			  z -= S(-1,t2,t4)*mu[t3]*cR2;
			  z *= -2.*eta[t1]*cR1;
			  break;
  case 3	: z  = S(+1,t1,t4)*S(-1,t2,t3)*cR1*cL2;
			  z -= mu[t1]*mu[t2]*eta[t3]*eta[t4]*cL1*cL2;
			  z -= mu[t3]*mu[t4]*eta[t1]*eta[t2]*cR1*cR2;
			  z *= -2.;
			  break;
  case 4	: z  = S(+1,t3,t1)*mu[t2]*cR1;
			  z -= S(+1,t3,t2)*mu[t1]*cL1;
			  z *= -2.*eta[t4]*cR2;
			  break;
  case 5	: z = ( 0., 0. );
			  break;
  case 6	: z  = mu[t1]*mu[t4]*eta[t2]*eta[t3]*cL1*cL2;
			  z += mu[t2]*mu[t3]*eta[t1]*eta[t4]*cR1*cR2;
			  z -= mu[t1]*mu[t3]*eta[t2]*eta[t4]*cL1*cR2;
			  z -= mu[t2]*mu[t4]*eta[t1]*eta[t3]*cR1*cL2;
			  z *= -2.;
  case 7	: z  = S(+1,t2,t4)*mu[t1]*cL1;
			  z -= S(+1,t1,t4)*mu[t2]*cR1;
			  z *= -2.*eta[t3]*cL2;
			  break;

  case 8	: z  = S(-1,t2,t4)*mu[t1]*cR1;
  			  z -= S(-1,t1,t4)*mu[t2]*cL1;
  			  z *= -2.*eta[t3]*cR2;
  			  break;
  case 9	: z  = mu[t1]*mu[t4]*eta[t2]*eta[t3]*cR1*cR2;
			  z += mu[t2]*mu[t3]*eta[t1]*eta[t4]*cL1*cL2;
			  z -= mu[t1]*mu[t3]*eta[t2]*eta[t4]*cR1*cL2;
			  z -= mu[t2]*mu[t4]*eta[t1]*eta[t3]*cL1*cR2;
			  z *= -2.;
  case 10	: z = ( 0., 0. );
			  break;
  case 11	: z  = S(-1,t3,t1)*mu[t2]*cL1;
			  z -= S(-1,t3,t2)*mu[t1]*cR1;
			  z *= -2.*eta[t4]*cL2;
			  break;
  case 12	: z  = S(-1,t1,t4)*S(+1,t2,t3)*cL1*cR2;
			  z -= mu[t1]*mu[t2]*eta[t3]*eta[t4]*cR1*cR2;
			  z -= mu[t3]*mu[t4]*eta[t1]*eta[t2]*cL1*cL2;
			  z *= -2.;
			  break;
  case 13	: z  = S(+1,t2,t3)*mu[t4]*cR2;
			  z -= S(+1,t2,t4)*mu[t3]*cL2;
			  z *= -2.*eta[t1]*cL1;
			  break;
  case 14	: z  = S(-1,t4,t1)*mu[t3]*cR2;
			  z -= S(-1,t3,t1)*mu[t4]*cL2;
			  z *= -2.*eta[t2]*cL1;
			  break;
  case 15	: z  = S(-1,t3,t1)*S(+1,t4,t2)*cL1*cL2;
			  z -= mu[t1]*mu[t2]*eta[t3]*eta[t4]*cR1*cL2;
			  z -= mu[t3]*mu[t4]*eta[t1]*eta[t2]*cL1*cR2;
			  z *= -2.;
			  break;
  }
  return z;
}  

Complex XYZFunc::Y( const int t1, const int t2, const int hel_comb, const Complex cR, const Complex cL )
{
  Complex y(0.,0.);
  switch( hel_comb ) {
	case 0:	y = cR*mu[t1]*eta[t2]+cL*mu[t2]*eta[t1];
			break;
	case 1:	y = cL*S(+1, t1, t2);
			break;
	case 2:	y = cR*S(-1, t1, t2);
			break;
	case 3:	y = cL*mu[t1]*eta[t2]+cR*mu[t2]*eta[t1];
			break;
  }
  return y;
}  

Complex XYZFunc::X( const int t1, const int t2, const int t3, 
	const int hel_comb, const Complex cR, const Complex cL )
{
  Complex x(0., 0.);
  switch( hel_comb ) {
	case 0:	x  = mu[t1]*mu[t3]*eta[t2]*eta[t2]*cL;
			x += mu[t2]*mu[t2]*eta[t1]*eta[t3]*cR;
			x += cR*S(+1,t1,t2)*S(-1,t2,t3); 
			break; 
	case 1:	x  = cL*mu[t1]*eta[t2]*S(+1,t2,t3);
			x += cR*mu[t3]*eta[t2]*S(+1,t1,t2);
			break; 
	case 2:	x  = cR*mu[t1]*eta[t2]*S(-1,t2,t3);
			x += cL*mu[t3]*eta[t2]*S(-1,t1,t2);
			break; 
	case 3:	x  = mu[t1]*mu[t3]*eta[t2]*eta[t2]*cR;
			x += mu[t2]*mu[t2]*eta[t1]*eta[t3]*cL;
			x += cL*S(+1,t1,t2)*S(-1,t2,t3); 
			break;  
  }
  return x;
}  
 
Complex XYZFunc::Z( 
	const int t1, const int l1, 
	const int t2, const int l2,
	const int t3, const int l3,
	const int t4, const int l4,
	const Complex cR1, const Complex cL1, const Complex cR2, const Complex cL2 ) 
{											// l=0 <=> -; l=1 <=> +; <---- helicity
  const int hel_comb = ((1-l1)<<3) + ((1-l2)<<2) + ((1-l3)<<1) + (1-l4);
  return Z(t1,t2,t3,t4,hel_comb,cR1,cL1,cR2,cL2);
}
 
Complex XYZFunc::Y( 
	const int t1, const int l1, 
	const int t2, const int l2,
	const int t3, const int l3,
	const Complex cR, const Complex cL ) 
{											// l=0 <=> -; l=1 <=> +; <---- helicity
  const int hel_comb = ((1-l1)<<2) + ((1-l2)<<1) + (1-l3);
  return( Y(t1,t2,hel_comb,cR,cL) );
}
 

Complex XYZFunc::X( 
	const int t1, const int l1, 
	const int t2, const int l2,
	const int t3, const int l3,
	const Complex cR, const Complex cL ) 
{											// l=0 <=> -; l=1 <=> +; <---- helicity
  const int hel_comb = ((1-l1)<<2) + ((1-l2)<<1) + (1-l3);
  return( X(t1,t2,t3,hel_comb,cR,cL) );
}
 
Complex XYZFunc::Q(short hel)
{
  Complex q(0.,0.);
  q  = Y(0,0,hel,Complex(1.,0.),Complex(1.,0.)); 
  q -= X(0,N,0,hel,Complex(1.,0.),Complex(-1.,0.));	
  return q/sqrt(p_flav[0].PSMass()-SQRT_05*mu[0]*eta[0])/2./mu[0];
}
