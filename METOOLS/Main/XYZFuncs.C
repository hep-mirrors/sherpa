#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"


#ifndef SQRT_05
#define SQRT_05 0.70710678118654757
#endif

#ifndef GAUGEK0
#define GAUGEK0 10
#endif

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

// constructor
XYZFunc::XYZFunc(const Vec4D_Vector& p, const Flavour_Vector& fl, 
                 bool anti, const vector<int>& indices) :
  m_k0n(GAUGEK0), m_anti(anti)
{
  m_anti=anti;
  if (indices.size()==0) {
    p_i=NULL;
    m_N=fl.size();
  }
  else {
    p_i=&indices.front();
    m_N = indices.size();
  }
  p_mom = new Vec4D[m_N];
  p_flav = new Flavour[m_N];
  for( int i=0; i<m_N; i++ ) {
    p_mom[i] = p_i?p[p_i[i]]:p[i];
    p_flav[i] = p_i?fl[p_i[i]]:fl[i];
  }
  CalcEtaMu();
}

XYZFunc::XYZFunc( int n, const Vec4D* p, const Flavour *fl,
                  bool anti, const int *indices ) :
  m_N(n), m_k0n(GAUGEK0), m_anti(anti), p_i(indices)
{
  m_anti=anti;
  m_N = n;
  p_mom = new Vec4D[m_N];
  p_flav = new Flavour[m_N];
  for( int i=0; i<m_N; i++ ) {
    p_mom[i] = p_i?p[p_i[i]]:p[i];
    p_flav[i] = p_i?fl[p_i[i]]:fl[i];
  }
  CalcEtaMu();
}

XYZFunc::XYZFunc(const Flavour_Vector& fl,
                 const vector<int>& indices ) :
  m_k0n(GAUGEK0)
{
  if (indices.size()==0) {
    p_i=NULL;
    m_N=fl.size();
  }
  else {
    p_i=&indices.front();
    m_N = indices.size();
  }

  p_mom  = new Vec4D[m_N];
  p_flav = new Flavour[m_N];
  for( int i=0; i<m_N; i++ ) p_flav[i] = p_i?fl[p_i[i]]:fl[i];
}

XYZFunc::~XYZFunc()
{
  delete [] p_mom;
  delete [] p_flav;
}

void XYZFunc::Prepare( const Vec4D_Vector& p, const bool anti )
{
  m_anti = anti;
  for( int i=0; i<m_N; i++ ) p_mom[i] = p_i?p[p_i[i]]:p[i];
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
    case 0 :
      _m_eta = csqrt( 2.*(pi[0]-(pi[1]+pi[3])*SQRT_05) );
      break;
    case 1  :
      _m_eta = csqrt( 2.*(pi[0]-(pi[2]+pi[3])*SQRT_05) );
      break;
    case 2  :
      _m_eta = csqrt( 2.*(pi[0]-(pi[1]+pi[2])*SQRT_05) );
      break;
    case 10 :
      _m_eta = csqrt( 2.*(pi[0]+pi[Spinor<double>::R3()]) );
      break;
    default :
        THROW(not_implemented, "");
    }
    m_eta.push_back( _m_eta );
    Complex help( p_flav[i].HadMass(), 0. );
    m_mu.push_back( help/m_eta[i] );
    if((p_flav[i].IsAnti() && m_anti==false) ||
       (!p_flav[i].IsAnti() && m_anti==true)) m_mu[i] *= -1.;
  }
}


// building blocks

Complex XYZFunc::S( const int s, const int i, const int j )
{
  Complex A = m_eta[j]/m_eta[i];
  Complex ret(0.,0.);
  Vec4D pi(p_mom[i]);
  Vec4D pj(p_mom[j]);
  switch( m_k0n ) {
  case 1 :
    ret = Complex( 1.*(double)s*pi[1], SQRT_05*(pi[2]-pi[3]) )*A;
    ret -= Complex( 1.*(double)s*pj[1], SQRT_05*(pj[2]-pj[3]) )/A;
    break;
  case 10:
    ret = Complex( 1.*(double)s*pi[Spinor<double>::R1()], -pi[Spinor<double>::R2()] )*A;
    ret -= Complex( 1.*(double)s*pj[Spinor<double>::R1()], -pj[Spinor<double>::R2()] )/A;
    break;
  default:
    THROW(not_implemented, "k0n choice not fully implemented yet.");
  }
  return ret;
}

Complex XYZFunc::S( const int s, const int i, const Vec4C p, Complex eta )
{
  Complex A = eta/m_eta[i];
  Complex ret(0.,0.);
  Vec4D pi = p_mom[i];
  Vec4C pj = p;
  Complex ci(0.0,1.0);
  switch( m_k0n ) {
  case 1 :
    ret  = (1.*(double)s*pi[1] + ci*SQRT_05*(pi[2]-pi[3]) )*A
          -(1.*(double)s*pj[1] + ci*SQRT_05*(pj[2]-pj[3]) )/A;
    break;
  case 10:
    ret  = (1.*(double)s*pi[Spinor<double>::R1()] - ci*pi[Spinor<double>::R2()] )*A
          -(1.*(double)s*pj[Spinor<double>::R1()] - ci*pj[Spinor<double>::R2()] )/A;
    break;
  default:
    THROW(not_implemented, "k0n choice not fully implemented yet.");
  }
  return ret;
}

Complex XYZFunc::S( const int s, const Vec4C p, Complex eta, const int j )
{
  Complex A = m_eta[j]/eta;
  Complex ret(0.,0.);
  Vec4C pi = p;
  Vec4D pj = p_mom[j];
  Complex ci(0.0,1.0);
  switch( m_k0n ) {
  case 1 :
    ret  = (1.*(double)s*pi[1] + ci*SQRT_05*(pi[2]-pi[3]) )*A
          -(1.*(double)s*pj[1] + ci*SQRT_05*(pj[2]-pj[3]) )/A;
    break;
  case 10:
    ret  = (1.*(double)s*pi[Spinor<double>::R1()] - ci*pi[Spinor<double>::R2()] )*A
          -(1.*(double)s*pj[Spinor<double>::R1()] - ci*pj[Spinor<double>::R2()] )/A;
    break;
  default:
    THROW(not_implemented, "k0n choice not fully implemented yet.");
  }
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
  case 5	: z  = Complex( 0., 0. );
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
  case 10	: z  = Complex( 0., 0. );
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

Complex XYZFunc::X( const int t1, const Vec4C p2, const int t3, 
		    const int hel_comb, const Complex cR, const Complex cL )
{
  Complex x(0., 0.);
  Complex eta2 (0.,0.);
  switch( m_k0n ) {
  case 1  : eta2 = sqrt( 2.*(p2[0]-(p2[2]+p2[3])*SQRT_05) );
    break;
  case 2  : eta2 = sqrt( 2.*(p2[0]-(p2[1]+p2[2])*SQRT_05) );
    break;
  case 10 : eta2 = sqrt( 2.*(p2[0]+p2[Spinor<double>::R3()]) );
    break;
  default : eta2 = sqrt( 2.*(p2[0]-(p2[1]+p2[3])*SQRT_05) );
  }
  Complex mu2 = sqrt(p2*p2)/eta2;
  
  switch( hel_comb ) {
  case 0:
    x  = m_mu[t1]*m_mu[t3]*eta2*eta2*cL;
    x += mu2*mu2*m_eta[t1]*m_eta[t3]*cR;
    x += cR*S(+1,t1,p2,eta2)*S(-1,p2,eta2,t3);
    break;
  case 1:
    x  = cL*m_mu[t1]*eta2*S(+1,p2,eta2,t3);
    x += cR*m_mu[t3]*eta2*S(+1,t1,p2,eta2);
    break;
  case 2:
    x  = cR*m_mu[t1]*eta2*S(-1,p2,eta2,t3);
    x += cL*m_mu[t3]*eta2*S(-1,t1,p2,eta2);
    break;
  case 3:
    x  = m_mu[t1]*m_mu[t3]*eta2*eta2*cR;
    x += mu2*mu2*m_eta[t1]*m_eta[t3]*cL;
    x += cL*S(-1,t1,p2,eta2)*S(+1,p2,eta2,t3);
    break;
  }
  return x;
}

Vec4C XYZFunc::L( const int t1, const int t3,
                         const int hel_comb, const Complex cR, const Complex cL )
{
  Vec4C l( Complex(0.,0.), Complex(0.,0.), Complex(0.,0.), Complex(0.,0.) );
  ATOOLS::Vec4D k0, k1;
  switch( m_k0n ) {
  case 0  :
    k0 = ATOOLS::Vec4D(1., SQRT_05, 0., SQRT_05);
    k1 = ATOOLS::Vec4D(0., 0., 1., 0.);
    break;
  case 1 :
    k0 = ATOOLS::Vec4D(1., 0., SQRT_05, SQRT_05);
    k1 = ATOOLS::Vec4D(0., 1., 0., 0.);
    break;
  case 2  :
    k0 = ATOOLS::Vec4D(1., SQRT_05, SQRT_05, 0.);
    k1 = ATOOLS::Vec4D(0., 0., 0., 1.);
    break;
  case 10:
    k0 = Spinor<double>::GetK0();
    k1 = Spinor<double>::GetK1();
    break;
  default :
      THROW(not_implemented, "");
  }
  
  Complex A, B;
  ATOOLS::Vec4D tmp_cross, tmp_cross2;
  switch( hel_comb ) {
    case 0: // ++
      A = 2.*cR/m_eta[t1]/m_eta[t3];
      tmp_cross = cross(k0,p_mom[t1],p_mom[t3]);
      for(int i=0; i<4; i++) {
        l[i] = A*( (k0*p_mom[t1])*p_mom[t3][i] + (k0*p_mom[t3])*p_mom[t1][i] - (p_mom[t1]*p_mom[t3])*k0[i] + Complex(0.0,1.0)*tmp_cross[i] );
        l[i] += 2.0*cL*m_mu[t1]*m_mu[t3]*k0[i];
      }
      break;
    case 1: // +-
      A = 2.0*cR*m_mu[t3]/m_eta[t1];
      B = 2.0*cL*m_mu[t1]/m_eta[t3];
      tmp_cross = cross( k1, k0, p_mom[t1]);
      tmp_cross2 = cross( k1, k0, p_mom[t3]);
      for(int i=0; i<4; i++) {
        l[i] = A*( -k0[i]*(k1*p_mom[t1]) + k1[i]*(k0*p_mom[t1]) + Complex(0.0,1.0)*tmp_cross[i] )
              +B*( -k1[i]*(k0*p_mom[t3]) + k0[i]*(k1*p_mom[t3]) - Complex(0.0,1.0)*tmp_cross2[i] );
      }
      break;
    case 2: // -+
      A = -2.0*cL*m_mu[t3]/m_eta[t1];
      B = -2.0*cR*m_mu[t1]/m_eta[t3];
      tmp_cross = cross( k1, k0, p_mom[t1]);
      tmp_cross2 = cross( k1, k0, p_mom[t3]);
      for(int i=0; i<4; i++) {
        l[i] = A*( -k0[i]*(k1*p_mom[t1]) + k1[i]*(k0*p_mom[t1]) - Complex(0.0,1.0)*tmp_cross[i] )
            +B*( -k1[i]*(k0*p_mom[t3]) + k0[i]*(k1*p_mom[t3]) + Complex(0.0,1.0)*tmp_cross2[i] );
      }
      break;
    case 3: // ++
      A = 2.*cL/m_eta[t1]/m_eta[t3];
      tmp_cross = cross(k0,p_mom[t1],p_mom[t3]);
      for(int i=0; i<4; i++) {
        l[i] = A*( (k0*p_mom[t1])*p_mom[t3][i] + (k0*p_mom[t3])*p_mom[t1][i] - (p_mom[t1]*p_mom[t3])*k0[i] - Complex(0.0,1.0)*tmp_cross[i] );
        l[i] += 2.0*cR*m_mu[t1]*m_mu[t3]*k0[i];
      }
      break;
  }
  return l;
}
 
Complex XYZFunc::Z( 
	const int t1, const int l1, 
	const int t2, const int l2,
	const int t3, const int l3,
	const int t4, const int l4,
	const Complex cR1, const Complex cL1, const Complex cR2, const Complex cL2 ) 
{
  const int hel_comb = ((1-l1)<<3) + ((1-l2)<<2) + ((1-l3)<<1) + (1-l4);
  return m_anti?
    conj(Z(t1,t2,t3,t4,hel_comb,conj(cR1),conj(cL1),conj(cR2),conj(cL2))):
    Z(t1,t2,t3,t4,hel_comb,cR1,cL1,cR2,cL2);
}
 
Complex XYZFunc::Y( 
	const int t1, const int l1,
	const int t2, const int l2,
	const Complex cR, const Complex cL ) 
{
  const int hel_comb = ((1-l1)<<1) + (1-l2);
  return m_anti?
    conj(Y(t2,t1,hel_comb,conj(cL),conj(cR))):
    Y(t2,t1,hel_comb,cR,cL);
}

Complex XYZFunc::X( 
                    const int t1, const int l1,
                    const Vec4C p2,
                    const int t3, const int l3,
                    const Complex cR, const Complex cL )
{
  const int hel_comb = ((1-l1)<<1) + (1-l3);
  return m_anti?
    conj(X(t3,p2,t1,hel_comb,conj(cR),conj(cL))):
    X(t3,p2,t1,hel_comb,cR,cL);
}

Complex XYZFunc::G(
                   const int t1, const int l1,
                   const Vec4C p2,
                   const int t3, const int l3 )
{
  Complex i(0.0,1.0);
  Complex cL(1.0,0.0);
  Complex cR(1.0,0.0);
  return i*((p2*(p_mom[t1]+p_mom[t3]))*Y(t1,l1,t3,l3,cL,cR)-
            (p_flav[t1].HadMass()+p_flav[t3].HadMass())*X(t1,l1,p2,t3,l3,cL,cR));
}
 
Vec4C XYZFunc::L(
	const int t1, const int l1,
	const int t2, const int l2,
	const Complex cR, const Complex cL ) 
{
  const int hel_comb = ((1-l1)<<1) + (1-l2);
  return m_anti?
    conj(L(t1,t2,hel_comb,conj(cR),conj(cL))):
    L(t1,t2,hel_comb,cR,cL);
}

int XYZFunc::MToL(int m)
{
  switch(m) {
  case -1:
    return 1;
  case 0:
    return 2;
  case 1:
    return 0;
  default:
    msg_Error()<<METHOD<<" index out of bounds: m="<<m<<endl;
    THROW(fatal_error, "Aborting.");
  }
  return -1;
}

#define CLEBSCH_GORDAN_SUM3(STRUCTURE, PTHREE, LTHREE, CONTRACTION) \
  Vec4C eps;\
  int m;\
  double CG;\
  switch(LTHREE) {\
  case 0:\
    m   = -1;\
    CG  = 1.0;\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    break;\
  case 1:\
    m   = -1;\
    CG  = sqrt(1.0/3.0);\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    m   = 0;\
    CG  = sqrt(2.0/3.0);\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    break;\
  case 2:\
    m   = 0;\
    CG  = sqrt(2.0/3.0);\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    m   = 1;\
    CG  = sqrt(1.0/3.0);\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    break;\
  case 3:\
    m   = 1;\
    CG  = 1.0;\
    eps = Polarization_Vector(p_mom[PTHREE])[MToL(m)];\
    if(!m_anti) eps=conj(eps);\
    result += CONTRACTION*CG*STRUCTURE;\
    break;\
  default:\
    msg_Error()<<METHOD<<" index out of bounds: l1="<<l1<<endl;\
    THROW(fatal_error, "Aborting.");\
  }


Vec4C XYZFunc::Y31(const int t1, const int l1,
                    const int t2, const int l2,
                    Complex cR, Complex cL )
{
  Vec4C result;
  CLEBSCH_GORDAN_SUM3(Y(t1,l1-m, t2,l2, cR, cL),t1,l1,eps)
  return result;
}


Vec4C XYZFunc::Y13(const int t1, const int l1,
                    const int t2, const int l2,
                    Complex cR, Complex cL )
{
  Vec4C result;

  CLEBSCH_GORDAN_SUM3(Y(t1,l1, t2,l2-m, cR, cL),t2,l2,eps)
  return result;
}


Vec4C XYZFunc::X31(const int t1, const int l1,
                    const ATOOLS::Vec4C p2,
                    const int t3, const int l3,
                    Complex cR, Complex cL )
{
  Vec4C result;
  CLEBSCH_GORDAN_SUM3(X(t1,l1-m, p2, t3,l3, cR, cL),t1,l1,eps)
  return result;
}


Vec4C XYZFunc::X13(const int t1, const int l1,
                    const ATOOLS::Vec4C p2,
                    const int t3, const int l3,
                    Complex cR, Complex cL )
{
  Vec4C result;
  CLEBSCH_GORDAN_SUM3(X(t1,l1, p2, t3,l3-m, cR, cL),t3,l3,eps)
  return result;
}


Vec4C XYZFunc::L31(const int t1, const int l1,
                    const ATOOLS::Vec4C p,
                    const int t2, const int l2,
                    Complex cR, Complex cL )
{
  Vec4C result;
  CLEBSCH_GORDAN_SUM3(L(t1,l1-m, t2,l2, cR, cL),t1,l1,p*eps)
  return result;
}


Vec4C XYZFunc::L13(const int t1, const int l1,
                    const ATOOLS::Vec4C p,
                    const int t2, const int l2,
                    Complex cR, Complex cL )
{
  Vec4C result;
  CLEBSCH_GORDAN_SUM3(L(t1,l1, t2,l2-m, cR, cL),t2,l2,p*eps)
  return result;
}
