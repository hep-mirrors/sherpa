#include "Vector.H"
#include "Flavour.H"
#include "XYZFuncs.H"

using namespace HADRONS;
using namespace ATOOLS;

// constructor
XYZFunc::XYZFunc( int nout, const Vec4D *p, const Flavour *fl, int k0_n, bool anti )
{
  m_anti=anti;
  m_N = nout+1;
  p_mom = new Vec4D[m_N];
  p_flav = new Flavour[m_N];
  for( int i=0; i<m_N; i++ ) {
    p_mom[i] = p[i];
    p_flav[i] = fl[i];
  }
  m_k0n = k0_n;
  CalcEtaMu();
}
 
void XYZFunc::CalcEtaMu() 
{
  Vec4D pi;
  m_mu.clear(); m_eta.clear();
  for( int i=0; i<m_N; i++ )
  {
    pi = p_mom[i];
    Complex _m_eta (0.,0.);
    switch( m_k0n ) {
      case 1  : _m_eta = csqrt( 2.*(pi[0]-(pi[2]+pi[3])*SQRT_05) );
                break;
      case 2  : _m_eta = csqrt( 2.*(pi[0]-(pi[1]+pi[2])*SQRT_05) );
                break;
      default : _m_eta = csqrt( 2.*(pi[0]-(pi[1]+pi[3])*SQRT_05) );
    }
    m_eta.push_back( _m_eta );
    Complex help( p_flav[i].PSMass(), 0. );
    m_mu.push_back( help/m_eta[i] );
    if(     p_flav[i].IsAnti() && m_anti==false
        || !p_flav[i].IsAnti() && m_anti==true  ) m_mu[i] *= -1.;
  }
}


// building blocks

Complex XYZFunc::S( const int s, const int i, const int j )
{
  Complex A = m_eta[j]/m_eta[i];
  Complex ret(0.,0.);
  Vec4D pi, pj;
  pi = p_mom[i];
  pj = p_mom[j];
  ret  = Complex( 1.*(double)s*pi[1], SQRT_05*(pi[2]-pi[3]) )*A;
  ret -= Complex( 1.*(double)s*pj[1], SQRT_05*(pj[2]-pj[3]) )/A;
  return ret;
}  

Complex XYZFunc::S( const int s, const int i, const Vec4D p, Complex eta )
{
  Complex A = eta/m_eta[i];
  Complex ret(0.,0.);
  Vec4D pi, pj;
  pi = p_mom[i];
  pj = p;
  ret  = Complex( 1.*(double)s*pi[1], SQRT_05*(pi[2]-pi[3]) )*A;
  ret -= Complex( 1.*(double)s*pj[1], SQRT_05*(pj[2]-pj[3]) )/A;
  return ret;
}  
 
Complex XYZFunc::S( const int s, const Vec4D p, Complex eta, const int j )
{
  Complex A = m_eta[j]/eta;
  Complex ret(0.,0.);
  Vec4D pi, pj;
  pi = p;
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
			  z -= m_mu[t1]*m_mu[t2]*m_eta[t3]*m_eta[t4]*cL1*cR2;
			  z -= m_mu[t3]*m_mu[t4]*m_eta[t1]*m_eta[t2]*cR1*cL2;
			  z *= -2.;
			  break;
  case 1	: z  = S(1,t4,t1)*m_mu[t3]*cL2;
			  z -= S(+1,t3,t1)*m_mu[t4]*cR2;
			  z *= -2.*m_eta[t2]*cR1;
			  break;
  case 2	: z  = S(-1,t2,t3)*m_mu[t4]*cL2;
			  z -= S(-1,t2,t4)*m_mu[t3]*cR2;
			  z *= -2.*m_eta[t1]*cR1;
			  break;
  case 3	: z  = S(+1,t1,t4)*S(-1,t2,t3)*cR1*cL2;
			  z -= m_mu[t1]*m_mu[t2]*m_eta[t3]*m_eta[t4]*cL1*cL2;
			  z -= m_mu[t3]*m_mu[t4]*m_eta[t1]*m_eta[t2]*cR1*cR2;
			  z *= -2.;
			  break;
  case 4	: z  = S(+1,t3,t1)*m_mu[t2]*cR1;
			  z -= S(+1,t3,t2)*m_mu[t1]*cL1;
			  z *= -2.*m_eta[t4]*cR2;
			  break;
  case 5	: z = ( 0., 0. );
			  break;
  case 6	: z  = m_mu[t1]*m_mu[t4]*m_eta[t2]*m_eta[t3]*cL1*cL2;
			  z += m_mu[t2]*m_mu[t3]*m_eta[t1]*m_eta[t4]*cR1*cR2;
			  z -= m_mu[t1]*m_mu[t3]*m_eta[t2]*m_eta[t4]*cL1*cR2;
			  z -= m_mu[t2]*m_mu[t4]*m_eta[t1]*m_eta[t3]*cR1*cL2;
			  z *= -2.;
  case 7	: z  = S(+1,t2,t4)*m_mu[t1]*cL1;
			  z -= S(+1,t1,t4)*m_mu[t2]*cR1;
			  z *= -2.*m_eta[t3]*cL2;
			  break;

  case 8	: z  = S(-1,t2,t4)*m_mu[t1]*cR1;
  			  z -= S(-1,t1,t4)*m_mu[t2]*cL1;
  			  z *= -2.*m_eta[t3]*cR2;
  			  break;
  case 9	: z  = m_mu[t1]*m_mu[t4]*m_eta[t2]*m_eta[t3]*cR1*cR2;
			  z += m_mu[t2]*m_mu[t3]*m_eta[t1]*m_eta[t4]*cL1*cL2;
			  z -= m_mu[t1]*m_mu[t3]*m_eta[t2]*m_eta[t4]*cR1*cL2;
			  z -= m_mu[t2]*m_mu[t4]*m_eta[t1]*m_eta[t3]*cL1*cR2;
			  z *= -2.;
  case 10	: z = ( 0., 0. );
			  break;
  case 11	: z  = S(-1,t3,t1)*m_mu[t2]*cL1;
			  z -= S(-1,t3,t2)*m_mu[t1]*cR1;
			  z *= -2.*m_eta[t4]*cL2;
			  break;
  case 12	: z  = S(-1,t1,t4)*S(+1,t2,t3)*cL1*cR2;
			  z -= m_mu[t1]*m_mu[t2]*m_eta[t3]*m_eta[t4]*cR1*cR2;
			  z -= m_mu[t3]*m_mu[t4]*m_eta[t1]*m_eta[t2]*cL1*cL2;
			  z *= -2.;
			  break;
  case 13	: z  = S(+1,t2,t3)*m_mu[t4]*cR2;
			  z -= S(+1,t2,t4)*m_mu[t3]*cL2;
			  z *= -2.*m_eta[t1]*cL1;
			  break;
  case 14	: z  = S(-1,t4,t1)*m_mu[t3]*cR2;
			  z -= S(-1,t3,t1)*m_mu[t4]*cL2;
			  z *= -2.*m_eta[t2]*cL1;
			  break;
  case 15	: z  = S(-1,t3,t1)*S(+1,t4,t2)*cL1*cL2;
			  z -= m_mu[t1]*m_mu[t2]*m_eta[t3]*m_eta[t4]*cR1*cL2;
			  z -= m_mu[t3]*m_mu[t4]*m_eta[t1]*m_eta[t2]*cL1*cR2;
			  z *= -2.;
			  break;
  }
  return z;
}  

Complex XYZFunc::Y( const int t1, const int t2, const int hel_comb, const Complex cR, const Complex cL )
{
  Complex y(0.,0.);
  switch( hel_comb ) {
	case 0:	y = cR*m_mu[t1]*m_eta[t2]+cL*m_mu[t2]*m_eta[t1];
			break;
	case 1:	y = cL*S(+1, t1, t2);
			break;
	case 2:	y = cR*S(-1, t1, t2);
			break;
	case 3:	y = cL*m_mu[t1]*m_eta[t2]+cR*m_mu[t2]*m_eta[t1];
			break;
  }
  return y;
}  

Complex XYZFunc::X( const int t1, const int t2, const int t3, 
	const int hel_comb, const Complex cR, const Complex cL )
{
  Complex x(0., 0.);
  switch( hel_comb ) {
	case 0:	x  = m_mu[t1]*m_mu[t3]*m_eta[t2]*m_eta[t2]*cL;
			x += m_mu[t2]*m_mu[t2]*m_eta[t1]*m_eta[t3]*cR;
			x += cR*S(+1,t1,t2)*S(-1,t2,t3); 
			break; 
	case 1:	x  = cL*m_mu[t1]*m_eta[t2]*S(+1,t2,t3);
			x += cR*m_mu[t3]*m_eta[t2]*S(+1,t1,t2);
			break; 
	case 2:	x  = cR*m_mu[t1]*m_eta[t2]*S(-1,t2,t3);
			x += cL*m_mu[t3]*m_eta[t2]*S(-1,t1,t2);
			break; 
	case 3:	x  = m_mu[t1]*m_mu[t3]*m_eta[t2]*m_eta[t2]*cR;
			x += m_mu[t2]*m_mu[t2]*m_eta[t1]*m_eta[t3]*cL;
			x += cL*S(-1,t1,t2)*S(+1,t2,t3); 
			break;  
  }
  return x;
}  
 
Complex XYZFunc::X( const int t1, const Vec4D p2, const int t3, 
	const int hel_comb, const Complex cR, const Complex cL )
{
  if( p2.IsZero() ) return Complex(0.,0.);
  Complex x(0., 0.);
  Complex eta2 (0.,0.);
  switch( m_k0n ) {
    case 1  : eta2 = csqrt( 2.*(p2[0]-(p2[2]+p2[3])*SQRT_05) );
              break;
    case 2  : eta2 = csqrt( 2.*(p2[0]-(p2[1]+p2[2])*SQRT_05) );
              break;
    default : eta2 = csqrt( 2.*(p2[0]-(p2[1]+p2[3])*SQRT_05) );
  }
  Complex mu2 = csqrt(p2.Abs2())/eta2;
  switch( hel_comb ) {
	case 0:	x  = m_mu[t1]*m_mu[t3]*eta2*eta2*cL;
			x += mu2*mu2*m_eta[t1]*m_eta[t3]*cR;
			x += cR*S(+1,t1,p2,eta2)*S(-1,p2,eta2,t3); 
			break; 
	case 1:	x  = cL*m_mu[t1]*eta2*S(+1,p2,eta2,t3);
			x += cR*m_mu[t3]*eta2*S(+1,t1,p2,eta2);
			break; 
	case 2:	x  = cR*m_mu[t1]*eta2*S(-1,p2,eta2,t3);
			x += cL*m_mu[t3]*eta2*S(-1,t1,p2,eta2);
			break; 
	case 3:	x  = m_mu[t1]*m_mu[t3]*eta2*eta2*cR;
			x += mu2*mu2*m_eta[t1]*m_eta[t3]*cL;
			x += cL*S(-1,t1,p2,eta2)*S(+1,p2,eta2,t3); 
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
  if(m_anti) {
    const int hel_comb = ((1-l2)<<3) + ((1-l1)<<2) + ((1-l4)<<1) + (1-l3);
    return Z(t2,t1,t4,t3,hel_comb,cR1,cL1,cR2,cL2);
  }
  else {
    const int hel_comb = ((1-l1)<<3) + ((1-l2)<<2) + ((1-l3)<<1) + (1-l4);
    return Z(t1,t2,t3,t4,hel_comb,cR1,cL1,cR2,cL2);
  }
}
 
Complex XYZFunc::Y( 
	const int t1, const int l1,
	const int t2, const int l2,
	const Complex cR, const Complex cL ) 
{											// l=0 <=> -; l=1 <=> +; <---- helicity
  if(m_anti) {
    const int hel_comb = ((1-l2)<<1) + (1-l1);
    return( Y(t2,t1,hel_comb,cR,cL) );
  }
  else {
    const int hel_comb = ((1-l1)<<1) + (1-l2);
    return( Y(t1,t2,hel_comb,cR,cL) );
  }
}
 

Complex XYZFunc::X( 
	const int t1, const int l1, 
	const int t2,
	const int t3, const int l3,
	const Complex cR, const Complex cL ) 
{											// l=0 <=> -; l=1 <=> +; <---- helicity
  if(m_anti) {
    const int hel_comb = ((1-l3)<<1) + (1-l1);
    return( X(t3,t2,t1,hel_comb,cR,cL) );
  }
  else {
    const int hel_comb = ((1-l1)<<1) + (1-l3);
    return( X(t1,t2,t3,hel_comb,cR,cL) );
  }
}

Complex XYZFunc::X( 
                    const int t1, const int l1,
                    const Vec4D p2,
                    const int t3, const int l3,
                    const Complex cR, const Complex cL )
{                                                                                       // l=0 <=> -; l=1 <=> +; <---- helicity
  if(m_anti) {
    const int hel_comb = ((1-l3)<<1) + (1-l1);
    return( X(t3,p2,t1,hel_comb,cR,cL) );
  }
  else {
    const int hel_comb = ((1-l1)<<1) + (1-l3);
    return( X(t1,p2,t3,hel_comb,cR,cL) );
  }
}
 
Complex XYZFunc::Q(short hel)
{
  Complex q(0.,0.);
  q  = Y(0,0,hel,Complex(1.,0.),Complex(1.,0.)); 
  q -= X(0,m_N,0,hel,Complex(1.,0.),Complex(-1.,0.));	
  return q/sqrt(p_flav[0].PSMass()-SQRT_05*m_mu[0]*m_eta[0])/2./m_mu[0];
}
