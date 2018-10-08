#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! D_0(p12,p22,p32,p43,s12,s23;m12,m22,m32,m42)
/*!

      \            /
    p2 \    m3    / p3         p1 + p2 + p3 + p4 = 0
        \ ______ /
         |      |              s12 = (p1+p2)^2 = (p3+p4)^2
      m2 |      | m4           s23 = (p2+p3)^2 = (p1+p4)^2
         |______|
        /   m1   \
    p1 /          \ p4
      /            \
*/

Complex METOOLS::K(const double& z, const Complex& sm12, const Complex& sm22) {
  // Analytic continuation via z -> z+ieps => creates -ieps as imaginary part
  // Complex rat = (z-sqr(sm12-sm22))/(4.*sm12*sm22);
  // Complex root = sqrt((rat-1.)/rat);
  // Complex invopr = 1./(1.+root);
  if (z == sqr(sm12-sm22)) {
    return -1.;
  }
  else {
    // return -pow(invopr,2.)/rat;
    return (1.-sqrt(1.-4.*sm12*sm22/(z-sqr(sm12-sm22))))/(1.+sqrt(1.-4.*sm12*sm22/(z-sqr(sm12-sm22))));
  }
}

METOOLS::DivArrC
METOOLS::Master_Box(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
                    const double&  s12, const double&  s23,
                    const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
                    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! *******************************************************************************************
  //! four massless internal lines
  if (IsZero(m12) && IsZero(m22) && IsZero(m32) && IsZero(m42)) {
    //! all legs on-shell
    //! Ellis, Zanderighi Box.1
    //! D_0(0,0,0,0;s12,s23;0,0,0,0) = 1/(s12s23)*(4/epsIR2 + 2/epsIR*ln(mu4/s12s23)
    //!                                            + ln2(mu2/-s12) + ln2(mu2/-s23)        s12 < 0
    //!                                            - ln2(-s12/-s23) -pi2)                 s23 < 0
    if (IsZero(p12) && IsZero(p22) && IsZero(p32) && IsZero(p42))
      return 1./(s12*s23)*DivArrC(0.,
                                  2.*log(mu2*mu2/(s12*s23)),
                                  4.,
                                  sqr(log(mu2/(-s12))) + sqr(log(mu2/(-s23)))
                                  - sqr(log((-s12)/(-s23))) - M_PI*M_PI,
                                  0.,0.);
    //! one leg off-shell
    //! Ellis, Zanderighi Box.2
    //! D_0(0,0,0,p2;s12,s23;0,0,0,0) = 1/(s12s23)*(2/epsIR2 + 2/epsIR*ln(-mu2p2/s12s23)
    //!                                             + ln2(mu2/-s12) + ln2(mu2/-s23) - ln2(mu2/-p2)
    //!                                             - 2DiLog(1-p2/s12) - 2DiLog(1-p2/s23)
    //!                                             - ln2(-s12/-s23) - pi2/3)             s12 < 0
    //!                                                                                   s23 < 0
    //!                                                                                    p2 < 0
    //! continue to physical region by s12 -> s12 + i0
    //!                                s23 -> s23 + i0
    //!                                 p2 ->  p2 + i0
    //! maybe use other method -> easier to continue analytically
    else if (IsZero(p22) && IsZero(p32) && IsZero(p42)) 
      return Master_Box(p22,p32,p42,p12,s23,s12,m22,m32,m42,m12,mu2);
    else if (IsZero(p12) && IsZero(p32) && IsZero(p42)) 
      return Master_Box(p32,p42,p12,p22,s12,s23,m32,m42,m12,m22,mu2);
    else if (IsZero(p12) && IsZero(p22) && IsZero(p42)) 
      return Master_Box(p42,p12,p22,p32,s23,s12,m42,m12,m22,m32,mu2);
    else if (IsZero(p12) && IsZero(p22) && IsZero(p32)) {
      return 1./(s12*s23)*DivArrC(0.,
                                  2.*CLog(-mu2*p42/(s12*s23),-DivProdSign(p42,1,s12,1,s23,1)),
                                  2.,
                                  sqr(CLog(mu2/(-s12),1)) + sqr(CLog(mu2/(-s23),1)) - sqr(CLog(mu2/(-p42),1))
                                  - 2.*Dilog(1.-p42/s12,-DivSign(p42,1,s12,1)) - 2.*Dilog(1.-p42/s23,-DivSign(p42,1,s23,1))
                                  - sqr(CLog((-s12)/(-s23),DivSign(-s12,-1,-s23,-1))) - M_PI*M_PI/3.,
                                  0.,0.);
    }
    //! two opposite legs off-shell
    //! Ellis, Zanderighi Box.3
    //! D_0(0,p22,0,p42;s12,s23;0,0,0,0) = 1/(s12s23-p22p42)*(2/epsIR*ln(p22p42/s12s23)
    //!                                                       + ln2(mu2/-s12) + ln2(mu2/-s23)
    //!                                                       - ln2(mu2/-p22) - ln2(mu2/-p42)
    //!                                                       - 2DiLog(1-p22/s12) - 2DiLog(1-p22/s23)
    //!                                                       - 2DiLog(1-p42/s12) - 2DiLog(1-p42/s23)
    //!                                                       + 2DiLog(1-p22p42/s12s23)
    //!                                                       - ln2(-s12/-s23))
    //!                                                                         s12,s23,p22,p42 < 0
    //! continue to physical region by sij -> sij + i0
    //!                                pi2 -> pi2 + i0
    //! maybe use other method -> easier to continue analytically
    else if (IsZero(p12) && IsZero(p32))
      return 1./(s12*s23-p22*p42)*DivArrC(0.,
                                          2.*log(p22*p42/(s12*s23)),
                                          0.,
                                          sqr(log(mu2/(-s12))) + sqr(log(mu2/(-s23)))
                                          - sqr(log(mu2/(-p22))) - sqr(log(mu2/(-p42)))
                                          - 2.*DiLog(1.-p22/s12) - 2.*DiLog(1.-p22/s23)
                                          - 2.*DiLog(1.-p42/s12) - 2.*DiLog(1.-p42/s23)
                                          + 2.*DiLog(1.-(p22*p42)/(s12*s23))
                                          - sqr(log((-s12)/(-s23))),
                                          0.,0.);
    //! two adjacent legs off-shell
    //! Ellis, Zanderighi Box.4
    //! D_0(0,0,p32,p42;s12,s23;0,0,0,0) = 1/s12s23*(1/epsIR2 + 1/epsIR*(ln(-s12mu2/p32p42)
    //!                                                                  + 2*ln(p32p42/s12s23))
    //!                                             + ln2(mu2/-s12) + ln2(mu2/-s23)
    //!                                             - ln2(mu2/-p32) - ln2(mu2/-p42)
    //!                                             + 1/2*ln2(-s12mu2/p32p42)
    //!                                             - 2DiLog(1-p32/s23) - 2DiLog(1-p42/s23)
    //!                                             - ln2(s12/s23))
    //!                                                                         s12,s23,p32,p42 < 0
    //! continue to physical region by sij -> sij + i0
    //!                                pi2 -> pi2 + i0
    else if (IsZero(p12) && IsZero(p22))
      return 1./(s12*s23)*DivArrC(0.,
                                  log((-s12*mu2)/(p32*p42)) + 2.*log(p32*p42/(s12*s23)),
                                  1.,
                                  sqr(log(mu2/(-s12))) + sqr(log(mu2/(-s23)))
                                  - sqr(log(mu2/(-p32))) - sqr(log(mu2/(-p42)))
                                  - 0.5*sqr(log((-s12*mu2)/(p32*p42)))
                                  - 2.*DiLog(1.-p32/s23) - 2.*DiLog(1.-p42/s23)
                                  - sqr(log((-s12)/(-s23))),
                                  0.,0.);
    //! three legs off-shell
    //! Ellis, Zanderighi Box.5
    //! D_0(0,p22,p32,p42;s12,s23,0,0,0,0) = 1/(s12s23-p22p42)*(2/epsIR*ln(sqrt(p22*p42)/-s23)
    //!                                                        + ln2(mu2/-s12) + ln2(mu2/-s23)
    //!                                                        - ln2(mu2/-p22) - ln2(mu2/-p32)
    //!                                                        - ln2(mu2/-p42)
    //!                                                        + 1/2ln2(-s12mu2/p22p32)
    //!                                                        + 1/2ln2(-s12mu2/p32p42)
    //!                                                        - 2DiLog(1-p22/s12) - 2DiLog(1-p42/s23)
    //!                                                        - ln2(-s12/-s23))
    //!                                                                         s12,s23,p32,p42 < 0
    //! continue to physical region by sij -> sij + i0
    //!                                pi2 -> pi2 + i0
    else if (IsZero(p12))
      return 1./(s12*s23-p22*p42)*DivArrC(0.,
                                          log(sqrt(p22*p42)/(-s23)),
                                          0.,
                                          sqr(log(mu2/(-s12))) - sqr(log(mu2/(-s23)))
                                          - sqr(log(mu2/(-p22))) - sqr(log(mu2/(-p32)))
                                          - sqr(log(mu2/(-p42)))
                                          + 0.5*sqr(log((-s12*mu2)/(p22*p32)))
                                          + 0.5*sqr(log((-s12*mu2)/(p32*p42)))
                                          - 2.*DiLog(1.-p22/s12) - 2.*DiLog(1.-p42/s23)
                                          - sqr(log((-s12)/(-s23))),
                                          0.,0.);
    //! four external legs off-shell
    //! Ellis, Zanderighi Box.6-10
    //! D_0(p12,p22,p32,p43;s12,s23;0,0,0,0) =
    else
      return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  //! *******************************************************************************************
  //! three massless internal lines
    //! Ellis, Zanderighi Box.11-13
  else if (IsZero(m12) && IsZero(m22) && IsZero(m32))
    return DivArrC(0.,0.,0.,0.,0.,0.);
  //! *******************************************************************************************
  //! two massless internal lines
    //! Ellis, Zanderighi Box.14-15
  else if (IsZero(m12) && IsZero(m22))
    return DivArrC(0.,0.,0.,0.,0.,0.);
  //! *******************************************************************************************
  //! one massless internal line
    //! Ellis, Zanderighi Box.16
  else if (IsZero(m22) && !IsZero(m12) && !IsZero(m32) && !IsZero(m42)) 
    return Master_Box(p22,p32,p42,p12,s23,s12,m22,m32,m42,m12,mu2);
  else if (IsZero(m32) && !IsZero(m12) && !IsZero(m22) && !IsZero(m42)) 
    return Master_Box(p32,p42,p12,p22,s12,s23,m32,m42,m12,m22,mu2);
  else if (IsZero(m42) && !IsZero(m12) && !IsZero(m22) && !IsZero(m32)) 
    return Master_Box(p42,p12,p22,p32,s23,s12,m42,m12,m22,m32,mu2);

  else if (IsZero(m12) && IsZero(p22) && IsEqual(p12,p42) && IsEqual(p12,m42) && !IsEqual(p12,p32) && IsEqual(m22,m42) && IsEqual(m32,m42)) {
    Complex m2 = 1./3.*(m22+m32+m42);
    Complex b23 = sqrt(1.-4.*m42/s23);
    Complex x23 = (b23-1.)/(b23+1.);
    Complex b3 = sqrt(1.-4.*m42/p32);
    Complex x3 = (b3-1.)/(b3+1.);
    return 1./((s12-m2)*s23*b23)*DivArrC(0.,CLog(x23,1),0.,
					 CLog(x23,1)*CLog(mu2/m2,1)-Dilog(sqr(x23),1)-2.*CLog(x23,1)*CLog(1.-sqr(x23),-1)-2.*CLog(x23,1)*CLog(1.-s12/m2,-1)-sqr(CLog(x3,1))+sqr(M_PI)/6.-2.*(Dilog(1.-x23*x3,-Prod2Sign(x23,1,x3,1))+CLog(1.-x23*x3,-Prod2Sign(x23,1,x3,1))*(CLog(x23*x3,Prod2Sign(x23,1,x3,1))-CLog(x3,1)-CLog(x23,1)) + Dilog(1.-x23/x3,-DivSign(x23,1,x3,1))+CLog(1.-x23/x3,-DivSign(x23,1,x3,1))*(CLog(x23/x3,DivSign(x23,1,x3,1))-CLog(1./x3,-1)-CLog(x23,1))),0.,0.);
  }
  else if (IsZero(m12)) {
    Complex sm22 = sqrt(m22);
    Complex sm32 = sqrt(m32);
    Complex sm42 = sqrt(m42);
    Complex x23 = -K(s23,sm22,sm42);
    Complex x2 = -K(p22,sm22,sm32);
    Complex x3 = -K(p32,sm32,sm42);
    int ieps = 1; // -K(s,m1,m2) comes with +ieps or imaginary part of respective result
    // full expression:
    if (IsEqual(s23,sqr(sm22-sm42))) {
      return 1./(2.*sm22*sm42*(s12-m32))*
	DivArrC(0.,1.,0.,
		2.*CLog(sm32*sqrt(mu2)/(m32-s12),DivSign(sm32*sqrt(mu2),1,(m32-s12),-1)) 
		- (1.+x2*x3)/(1.-x2*x3)*(CLog(x2,ieps)+CLog(x3,ieps))
		- (x3+x2)/(x3-x2)*(CLog(x2,ieps)-CLog(x3,ieps)) - 2.,
		0.,0.);
    }
    else {
      return x23/(sm22*sm42*(s12-m32)*(1.-sqr(x23)))*
	DivArrC(0.,-CLog(x23,ieps),0.,
		-2.*CLog(x23,ieps)*CLog(sm32*sqrt(mu2)/(m32-s12),DivSign(sm32*sqrt(mu2),1,(m32-s12),-1)) 
		+ sqr(CLog(x2,ieps)) + sqr(CLog(x3,ieps))
		-Dilog(1.-sqr(x23),-ieps) + Dilog(1.-x23*x2*x3,-Prod3Sign(x23,ieps,x2,ieps,x3,ieps))
		+ Dilog(1.-x23/(x2*x3),-DivProdSign(x23,ieps,x2,ieps,x3,ieps))
		+ Dilog(1.-x23*x2/x3,-ProdDivSign(x23,ieps,x2,ieps,x3,ieps)) 
		+ Dilog(1.-x23*x3/x2,-ProdDivSign(x23,ieps,x3,ieps,x2,ieps)),
		0.,0.);
    }
  }
  //! *******************************************************************************************
  //! no massless internal lines
  else {
    //! all same mass
    if (IsEqual(m12,m22) && IsEqual(m22,m32) && IsEqual(m32,m42))
      //! D_0(0,0,0,0;s12,s23;m2,m2,m2,m2) = ?
      if (IsZero(p12) && IsZero(p22) && IsZero(p32) && IsZero(p42))
        return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}
