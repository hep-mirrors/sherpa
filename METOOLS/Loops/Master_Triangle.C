#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! C_0(p12,p22,p32;m12,m22,m32)
/*! all momenta incoming

        p2
       -----|\              p1 + p2 + p3 = 0
            | \ m3
            |  \            s12 = p3^2 = (-p1-p2)^2 = (p1+p2)^2
            |   \   p3
         m2 |    |-----
            |   /
            |  /            s12 > 0    s-channel
            | / m1          s12 < 0    t-channel
       -----|/
        p1
*/
METOOLS::DivArrC
METOOLS::Master_Triangle(const double&  p12, const double&  p22, const double&  p32,
                         const Complex& m12, const Complex& m22, const Complex& m32,
                         double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! three massless internal lines
  if (IsZero(m12) && IsZero(m22) && IsZero(m32)) {
    //! C_0(0,0,0;0,0,0) = ?
    if (IsZero(p12) && IsZero(p22) && IsZero(p32))
      return DivArrC(0.,0.,0.,0.,0.,0.); // wrong?
    //! C_0(0,0,s;0,0,0) = 1/s*(1/epsIR2 + 1/epsIR*ln(mu2/-s) + 1/2*ln2(mu2/-s))         s < 0
    //! no s-channel expression found
    else if (IsZero(p12) && IsZero(p22))
      return (1./p32)
              *DivArrC(0.,log(mu2/(-p32)),1.,0.5*sqr(log(mu2/(-p32))),0.,0.);
    //! C_0(0,p2,s;0,0,0) = 1/(p2-s)*(1/epsIR*ln(s/p2) + 1/2*(ln2(mu2/p2)-ln2(mu2/s)))
    //! no s-channel expression found
    else if (IsZero(p12)) {
      if (!IsZero(p22-p32))
        return 1./(p22-p32)*DivArrC(0.,log(p32/p22),0.,0.5*(sqr(log(mu2/(-p22)))
                                                           -sqr(log(mu2/(-p32)))),0.,0.);
      else
        return -1./p32*DivArrC(0.,1.,0.,log(mu2/(-p32)),0.,0.);    // needs to be checked
    }
  }
  //! two massless internal lines
  else if (IsZero(m12) && IsZero(m22)) {
    Complex m2(m32);
    //! C_0(0,0,0;0,0,m2) = 1/epsIR + ln(mu^2/m^2) + 1
    if (IsZero(p12) && IsZero(p22) && IsZero(p32))
      return DivArrC(0.,1.,0.,log(mu2/m32)+1.,0.,0.);
    else if (IsEqual(p22,m2) && IsEqual(p32,m2)) {
      //! C_0(0,m2,m2;0,0,m2) = -1/2m2*(1/epsIR + ln(mu2/m2) - 1)
      if (IsZero(p12))
        return -0.5/m2*DivArrC(0.,1.,0.,log(mu2/m2)-1.,0.,0.);
      //! C_0(s,m2,m2;0,0,m2) = 1/bs(2pi2/3 + 2DiLog(-lm/lp) + 1/2ln2(-lm/lp)              s < 0
      //! s-channel expression not found
      else {
        Complex b  = sqrt(1.-4.*m2/p12);
        Complex lp = 0.5*(1.+b);
        Complex lm = 0.5*(1.-b);
        return 1./(p12*b)*DivArrC(0.,0.,0.,2.*M_PI*M_PI/3.+2.*DiLog(-lm/lp)+0.5*sqr(log(-lm/lp)),0.,0.);
      }
    }
    //! C_0(0,m2,p2;0,0,m2) = -1/(m2-p2)*(1/epsIR2 + 1/epsIR(ln(m2/(m2-p2))+1/2ln(mu2/m2))
    //!                                   + ln(mu2/m2)*ln(m2/(m2-p2)) + 1/2 ln2(m2/(m2-p2))
    //!                                   + pi2/12 - DiLog(-p2/(m2-p2)))                    p2 < 0
    //! no expression for s-channel found yet
    else if (IsZero(p12) && (IsEqual(p22,m2) || IsEqual(p32,m2))) {
      Complex p2 = p22;
      if (IsEqual(p2,m2)) p2 = p32;
//       if (p2 < 0.)
        return 1./(p2-m2)*DivArrC(0.,
                                  log(m2/(m2-p2))+0.5*log(mu2/m2),
                                  0.5,
                                  log(mu2/m2)*log(m2/(m2-p2))+0.5*sqr(log(m2/(m2-p2)))
                                      +M_PI*M_PI/12.-DiLog(-p2/(m2-p2)),
                                  0.,0.);
//       else
//         return DivArrC(0.,0.,0.,0.,0.); // no expression found
    }
    //! C_0(0,p22,p32;0,0,m2) = 1/(p22-p32)*(1/epsIR*ln((m2-p32)/(m2-p22))
    //!                                      + 1/2*ln2(mu2/-p22) - 1/2*ln2(mu2/-p32)
    //!                                      + ln((m2-p22)/-p22)*ln(-p22(m2-p22)/(mu2*p22))
    //!                                      - ln((m2-p32)/-p32)*ln(-p32(m2-p32)/(mu2*p32))
    //!                                      - DiLog(m2/p22) + DiLog(m2/p32))
    else if (IsZero(p12))
      return 1./(p22-p32)*DivArrC(0.,
                                  log((m2-p32)/(m2-p22)),
                                  0.,
                                  0.5*sqr(log(mu2/(-p22)))-0.5*sqr(log(mu2/(-p32)))
                                  +log((m2-p22)/(-p22))*log(-p22*(m2-p22)/(mu2*m2))
                                  -log((m2-p32)/(-p32))*log(-p32*(m2-p32)/(mu2*m2))
                                  -DiLog(m2/p22)+DiLog(m2/p32),
                                  0.,0.);
  }
  //! one massless internal line
  else if (IsZero(m12)) {
    if (IsEqual(p32,m32) && IsEqual(m22,m32) && IsZero(p22)) {
      Complex m2 = 0.5*(m22+m32);
      //! C_0(m2,0,m2;0,m2,m2) = 1/2m2*(1/epsIR + ln(mu2/m2) - 2 + 2*ln 2)
      if (IsEqual(p12,m22))
        return 0.5/m2*DivArrC(0.,1.,0.,log(mu2/m2)-2.+2.*log(2.),0.,0.);
      //! C_0(s,0,m2;0,m2,m2) = -1/(m2-s)*(pi2/6 - DiLog(s/m2))
      else
        return -1./(m2-p12)*DivArrC(0.,0.,0.,M_PI*M_PI/6. - DiLog(p12/m2),0.,0.);
    }
    //! C_0(m2,s,m2;0,m2,m2) = 1/sb*(1/epsIR*ln x + ln x *ln(mu2/m2) + ln2((b-1)/2b) 
    //!                                           - 1/2*ln2 x - 2*DiLog((b-1)/2b) - pi2/6)     s < 0
    //! no s-channel expression found
    else if (IsEqual(p12,m22) && IsEqual(p32,m32) && IsEqual(m22,m32)) {
      Complex m2 = 0.5*(m22+m32);
      Complex b  = sqrt(1.-4.*m2/p22);
      Complex x  = (b-1.)/(b+1.);
      return 1./(p22*b)*DivArrC(0.,
                                log(x),
                                0.,
                                log(x)*log(mu2/m2)+sqr(log((b-1.)/(2.*b)))
                                -0.5*sqr(log(x))-2.*DiLog((b-1.)/(2.*b))-M_PI*M_PI/6.,
                                0.,0.);
    }
    //! C_0(m22,s,m32;0,m22,m32) = 1/2sb*(1/epsIR*ln(xp*xm) - ln(mu2/bs)ln(xp*xm)
    //!                                                     + 1/2ln2(-gp) - 1/2ln2(1-gp)
    //!                                                     - 1/2ln2(gm) + 1/2ln2(gm-1)
    //!                                                     - DiLog((1-gm)/b) - DiLog(gp/b)
    //!                                                     + DiLog((gp-1)/b) + DiLog(-gm/b))  s < 0
    //! no s-channel expression found
    else if (IsEqual(p12,m22) && IsEqual(p32,m32)) {
      Complex b  = sqrt(sqr(-p22+m22+m32)-4.*m22*m32)/p22;
      Complex g0 = (m32-m22+p22)/p22;
      Complex gp = 0.5*(g0+b);
      Complex gm = 0.5*(g0-b);
      Complex xp = (gp-1.)/gp;
      Complex xm = -gm/(1.-gm);
      return 0.5/(b*p22)*DivArrC(0.,
                                 log(xp*xm),
                                 0.,
                                 log(mu2/(b*p22))*log(xp*xm)
                                 + 0.5*sqr(log(-gp)) - 0.5*sqr(log(1.-gp))
                                 - 0.5*sqr(log(gm)) + 0.5*sqr(log(gm-1.))
                                 - DiLog((1.-gm)/b) - DiLog(gp/b)
                                 + DiLog((gp-1.)/b) + DiLog(-gm/b),
                                 0.,0.);
    }
  }
  //! no massless internal lines
  else {
    if (IsEqual(m12,m22) && IsEqual(m22,m32)) {
      Complex m2 = (m12+m22+m32)/3.;
      //! C_0(0,0,0;m2,m2,m2) = -1/2m2
      if (IsZero(p12) && IsZero(p22) && IsZero(p32))
        return -0.5/m2*DivArrC(0.,0.,0.,1.,0.,0.);
      //! C_0(0,0,s;m2,m2,m2) = 1/2s*ln2(-(1+b)/(1-b))                                        s < 0
      //! no s-channel expression found
      else if (IsZero(p12) && IsZero(p22)) {
        Complex b  = sqrt(1.-4.*m2/p32);
        return 0.5/p32*DivArrC(0.,0.,0.,sqr(log(-(1.+b)/(1.-b))),0.,0.);
      }
    }
    //! general finite box not implemented yet
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}


