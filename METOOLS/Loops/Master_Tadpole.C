#include "METOOLS/Loops/Master_Integrals.H"
#include "ATOOLS/Math/MathTools.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! A_0(m2)
/*!
            -------
          /        \    q = 0
       m |         |----------
         \        /
          -------
*/
METOOLS::DivArrC
METOOLS::Master_Tadpole(const Complex& m2, double mu2=0.) {

  DivArrC da(1.,1.,1.,1.,1.,1.);
  da*=Complex(0.,2.);

  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
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

