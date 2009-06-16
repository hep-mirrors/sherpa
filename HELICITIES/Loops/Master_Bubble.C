#include "HELICITIES/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace HELICITIES;

//! B_0(s;m02,m12)
/*!
            -------
      q   /   m1   \   q
    -----|         |------    s = q^2
         \   m2   /
          -------
*/

HELICITIES::DivArrC
HELICITIES::Master_Bubble(const double& s, 
                          const Complex& m02, const Complex& m12,
                          const double& mu2=0.) {
//   if (mu2 == 0)   mu2 = GLOBAL_RENORMALIZATION_SCALE;
  Complex m2 = 0.5*(m02+m12);
  //! two massless internal lines
  if (IsZero(m02) && IsZero(m12)) {
    //! B_0(0;0,0) = 1/epsUV - 1/epsIR
    if (IsZero(s))
      return DivArrC(1.,-1.,0.,0.,0.,0.);
    //! B_0(s;0,0) = 1/epsUV + ln(mu^2/s) + 2
    else if (IsZero(m02) && IsZero(m12))
      return DivArrC(1.,0.,0.,log(mu2/s)+2.,0.,0.);
  }
  //! one massless internal line
  else if (IsZero (m02) || IsZero(m12)) {
    //! B_0(0;0,m^2) = 1/epsUV + ln(mu^2/m^2) + 1
    if (IsZero(s))
      return DivArrC(1.,0.,0.,log(mu2/(2.*m2))+1.,0.,0.);
    //! B_0(s;0,m^2)
    else {
      //! a) B_0(m^2;0,m^2) = 1/epsUV + ln(mu^2/m^2) + 2
      if (IsEqual(s,m2))
        return DivArrC(1.,0.,0.,log(mu2/(2.*m2))+2.,0.,0.);
      //! b) B_0(s;0,m^2) = 1/epsUV + ln (mu^2/m^2) + (s+m^2)/s*ln((s+m^2)/s) + 2
      else
        return DivArrC(1.,0.,0.,log(mu2/m2) + 2. + (s+m2)/s*log(m2/(s+m2)),0.,0.);
    }
  }
  //! no massless internal line
  else {
    //! B_0(0;m^2,m^2) = 1/epsUV + ln(mu^2/m^2)
    if ((IsZero(m02-m12)) && (IsZero(s)))
      return DivArrC(1.,0.,0.,log(mu2/m2),0.,0.);
    //! B_0(s;m^2,m^2) = 1/epsUV + ln(mu^2/m^2) - beta*ln(-l+/l-) + 2
    else if (IsZero(m02-m12)) {
      Complex beta = sqrt(1.+4.*(m2/s));
      Complex lp   = 0.5*(1.+beta);
      Complex lm   = 0.5*(1.-beta);
      return DivArrC(1.,0.,0.,log(mu2/m2) + 2. - beta*log(-lp/lm),0.,0.);
    }
    //! B_0(0;m_1^2,m_2^2) = 1/epsUV + m02/(m02-m01)*ln(mu^2/m02) - m12/(m02-m12)*ln(mu2/m12) + 1
    else if (IsZero(s))
      return DivArrC(1.,0.,0.,m02/(m02-m12)*log(mu2/m02) - m12/(m02-m12)*log(mu2/m12) + 1.,0.,0.);
    //! B_0(s;m_1^2,m_2^2) = 1/epsUV + ln(mu^2/s) + sum(li*ln((li-1)/li) - ln(li-1)) + 2
    else {
      Complex l1 = (-s-m12+m02+sqrt(csqr(-s-m12+m02)+4.*s*m02))
                    /(-2.*s);
      Complex l2 = (-s-m12+m02-sqrt(csqr(-s-m12+m02)+4.*s*m02))
                    /(-2.*s);
      return DivArrC(1.,0.,0.,log(mu2/-s) + l1*log((l1-1.)/l1) + l2*log((l2-1.)/l2)
                                                               - log(l1-1.) - log(l2-1.) + 2.,0.,0.);
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}
