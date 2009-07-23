#include "HELICITIES/Loops/Master_Integrals.H"
#include "HELICITIES/Loops/PV_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace HELICITIES;

//! A_2(m2) in Passarino-Veltman reduction
/*!
            -------
          /        \    q = 0
       m |         |----------
         \        /
          -------
*/
HELICITIES::DivArrC
HELICITIES::PV_Tadpole_2(const Complex& m2, double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! massless internal line
  //! A_2(0) = 0
  if (IsZero(m2)) {
    return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  //! massive internal line
  //! A_2(m^2) = m^4(1/epsUV + ln(mu^2/m^2) + 3/2)
  else {
    return DivArrC(m2*m2,0.,0.,m2*m2*(log(mu2/m2) + 3./2.),0.,0.);
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

