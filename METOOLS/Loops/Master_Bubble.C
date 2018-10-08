#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! B_0(s;m02,m12)
/*!
            -------
      q   /   m1   \   q
    -----|         |------    s = q^2
         \   m2   /
          -------
*/

METOOLS::DivArrC
METOOLS::Master_Bubble(const double& s,
                       const Complex& m02, const Complex& m12,
                       double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  Complex m2 = 0.5*(m02+m12);
  //! two massless internal lines
  if (IsZero(m02) && IsZero(m12)) {
    //! B_0(0;0,0) = 1/epsUV - 1/epsIR
    if (IsZero(s) || s < pow(10.,-6.))
      return DivArrC(1.,-1.,0.,0.,0.,0.);
    //! B_0(s;0,0) = 1/epsUV + ln(mu^2/s) + 2
    else
      return DivArrC(1.,0.,0.,CLog(mu2/-s,1)+2.,0.,0.);
  }
  //! one massless internal line
  else if (IsZero (m02) || IsZero(m12)) {
    //! B_0(0;0,m^2) = 1/epsUV + ln(mu^2/m^2) + 1
    if (IsZero(s))
      return DivArrC(1.,0.,0.,CLog(mu2/(2.*m2),1)+1.,0.,0.);
    //! B_0(s;0,m^2)
    else {
      //! a) B_0(m^2;0,m^2) = 1/epsUV + ln(mu^2/m^2) + 2
      if (IsEqual(s,2.*m2))
        return DivArrC(1.,0.,0.,CLog(mu2/(2.*m2),1)+2.,0.,0.);
      //! b) B_0(s;0,m^2) = 1/epsUV + ln (mu^2/m^2) - (s+m^2)/s*ln((s+m^2)/m^2) + 2
      else {
	return DivArrC(1.,0.,0.,CLog(mu2/(2.*m2),1) + 2. + (2.*m2-s)/s*CLog((2.*m2-s)/(2.*m2),-1),0.,0.);}
    }
  }
  //! no massless internal line
  else {
    //! B_0(0;m^2,m^2) = 1/epsUV + ln(mu^2/m^2)
    if ((IsZero(m02-m12)) && (IsZero(s) || s < pow(10.,-6.)))
      return DivArrC(1.,0.,0.,log(mu2/m2),0.,0.);
    //! B_0(s;m^2,m^2) = 1/epsUV + ln(mu^2/m^2) - beta*ln(-l+/l-) + 2
    else if (IsZero(m02-m12)) {
      Complex beta = sqrt(1.-4.*(m2/s));
      Complex lp   = 0.5*(1.+beta);
      Complex lm   = 0.5*(1.-beta);
      Complex r = (-s+2.*m2+sqrt(sqr(s)-4.*m2*s))/(2.*m2);
      int ieps = ReSign(r-1./r);
      if (r == 0. && s > abs(m2)) return Master_Bubble(s,0.,0.,mu2);
      if (r == 0.) msg_Out() << "B_0(s,m2,m2,mu2) nan, s: " << s << " , m2: " << m2 << " , mu2: " << mu2 << "\n";
      // return DivArrC(1.,0.,0.,log(mu2/m2) + 2. - beta*log(-lp/lm),0.,0.);
      return DivArrC(1.,0.,0.,2.-CLog(m2/mu2,-1)-m2/s*(1./r-r)*CLog(r,ieps),0.,0.);
    }
    //! B_0(0;m_1^2,m_2^2) = 1/epsUV + m02/(m02-m01)*ln(mu^2/m02) - m12/(m02-m12)*ln(mu2/m12) + 1
    else if (IsZero(s) || s < pow(10.,-6.))
      return DivArrC(1.,0.,0.,m02/(m02-m12)*CLog(mu2/m02,1) - m12/(m02-m12)*CLog(mu2/m12,1) + 1.,0.,0.);
    //! B_0(s;m_1^2,m_2^2) = 1/epsUV + ln(mu^2/s) + sum(li*ln((li-1)/li) - ln(li-1)) + 2
    else {
      Complex r = (m02+m12-s+sqrt(csqr(m02+m12-s)-4.*m02*m12))/(2.*sqrt(m02)*sqrt(m12));
      int ieps = ReSign(r-1./r);
      // Denner eq. 4.23
      return DivArrC(1.,0.,0.,2. - CLog(sqrt(m02)*sqrt(m12)/mu2,-1) 
		     + (m02-m12)/s*CLog(sqrt(m12)/sqrt(m02),(std::real(m12)>std::real(m02))?-1:1)
		     - sqrt(m02)*sqrt(m12)/s*(1./r-r)*CLog(r,ieps),0.,0.);
    }
  }
  msg_Out() << "No bubble result found!\n"
	    << "p2: " << s << "\n"
	    << "(m12, m22): (" << m02 << ", " << m12 << ")\n";
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

// Derivative of B_0. Mostly based on direct calculation or Bardin:1999ak
METOOLS::DivArrC
METOOLS::Master_Bubble_Prime(const double& s,
			     const Complex& m02, const Complex& m12,
			     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (IsZero(s)) {
    if (IsZero(m02) && IsZero(m12)) {
      //msg_Out()<<"Call in wrong situation\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    if (IsZero(m12) || IsZero(m02)) {
      Complex m2 = m02+m12;
      return 0.5/m2*DivArrC(0.,0.,0.,1.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      Complex m2 = 0.5*(m02+m12);
      return DivArrC(0.,0.,0.,1./(6.*m2),0.,0.);
    }
    else {
      return 1./sqr(m02-m12)*DivArrC(0.,0.,0.,(m02+m12)/2.-m02*m12/(m02+m12)*CLog(m02/m12,(std::real(m02)>std::real(m12))?-1:1),0.,0.);
    }
  }
  if (IsZero(m02) && IsZero(m12)) {
    return 1./s*DivArrC(0.,0.,0.,-1.,0.,0.);
  }
  if (IsZero(m12) || IsZero(m02)) {
    Complex m2 = m02+m12;
    if (IsEqual(s,m2) || abs(m2/s) > pow(10.,6.)) {
      return 0.5/m2*DivArrC(0.,-1.,0.,-2.-CLog(mu2/m2,1),0.,0.);
    }
    return DivArrC(0.,0.,0.,-1./s+m2/sqr(s)*CLog(m2/(s-m2),-1),0.,0.);
  }
  if (IsEqual(m02,m12)) {
    Complex m2 = 0.5*(m02+m12);
    Complex beta = sqrt(1.-4.*(m2/s));
    return DivArrC(0.,0.,0.,-1./s+2.*m2*beta*CLog((2.*m2-s-s*beta)/(2.*m2),-1)/(s*(4.*m2-s)),0.,0.);
  }
  else {
    //! B_0_prime(s,m02,m12) = - (m02 - m12)/(2*s^2)*ln(m12/m02) + m0*m1/s^2*(1/r-r)*ln(r) - 1/s*(1+(r^2+1)/(r^2-1)*log(r)
    //! r = - 1/(2*m0*m1)*(sqrt((m02+m12-s)^2-4*m02*m12)-m02-m12+s)  // No issues here with analytic continuation since r does not cross negative real axis even for complex masses
    // Complex r = -0.5/(sqrt(m02)*sqrt(m12))*(sqrt(sqr(m02+m12-s)-4.*m02*m12)-m02-m12+s);
    // return DivArrC(0.,0.,0.,-0.5*(m02-m12)/sqr(s)*log(m12/m02) + sqrt(m02)*sqrt(m12)/sqr(s)*(1./r-r)*log(r) - 1./s*(1.+(sqr(r)+1.)/(sqr(r)-1.)*log(r)),0.,0.);
    if (IsEqual(s,m02) && abs(m12/m02) > pow(10.,6.)) return Master_Bubble_Prime(s,0.,m12,mu2);
    Complex l1 = (s-m12+m02+sqrt(csqr(s-m12+m02)-4.*s*m02))/(2.*s);
    Complex l2 = (s-m12+m02-sqrt(csqr(s-m12+m02)-4.*s*m02))/(2.*s);
    Complex dl1 = -(sqrt(csqr(s-m12+m02)-4.*s*m02)+s*(m02+m12-s)/sqrt(csqr(s-m12+m02)-4.*s*m02)+m02-m12)/(2.*sqr(s));
    Complex dl2 = (sqrt(csqr(s-m12+m02)-4.*s*m02)+s*(m02+m12-s)/sqrt(csqr(s-m12+m02)-4.*s*m02)-m02+m12)/(2.*sqr(s));
    // if (l1 == 0. || l2 == 0.) return Master_Bubble_Prime(s,0.,m02+m12,mu2);
    return DivArrC(0.,0.,0.,-1./s+dl1*log((l1-1.)/l1)+dl2*log((l2-1.)/l2),0.,0.);
  }
  msg_Out() << "No bubble derivative found!\n"
	    << "p2: " << s << "\n"
	    << "(m12, m22): (" << m02 << ", " << m12 << ")\n";
  return DivArrC(0.,0.,0.,0.,0.,0.);
}
