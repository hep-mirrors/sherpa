#include "HELICITIES/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace HELICITIES;

//! A_0(m2)
HELICITIES::DivArrC
HELICITIES::Master_Tadpole(const Complex& m2, const double& mu2=0.) {
//   if (mu2 == 0)   mu2 = GLOBAL_RENORMALIZATION_SCALE;
  //! massless internal line
  //! A_0(0) = 0
  if (IsZero(m2)) {
    return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  //! massive internal line
  //! A_0(m^2) = m^2(1/epsUV + ln(mu^2/m^2) + 1)
  else {
    return DivArrC(m2,0.,0.,m2*(log(mu2/m2) + 1.),0.,0.);
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

