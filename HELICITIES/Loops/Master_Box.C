#include "HELICITIES/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace HELICITIES;

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
HELICITIES::DivArrC
HELICITIES::Master_Box(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
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
    else if (IsZero(p12) && IsZero(p22) && IsZero(p32))
      return 1./(s12*s23)*DivArrC(0.,
                                  2.*log(-mu2*p42/(s12*s23)),
                                  2.,
                                  sqr(log(mu2/(-s12))) + sqr(log(mu2/(-s23))) - sqr(log(mu2/(-p42)))
                                  - 2.*DiLog(1.-p42/s12) - 2.*DiLog(1.-p42/s23)
                                  - sqr(log((-s12)/(-s23))) - M_PI*M_PI/3.,
                                  0.,0.);
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
  else if (IsZero(m12))
    return DivArrC(0.,0.,0.,0.,0.,0.);
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
