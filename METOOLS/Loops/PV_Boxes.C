#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! D_i(p12,p22,p32,p43,s12,s23;m12,m22,m32,m42)
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

METOOLS::DivArrC
METOOLS::PV_Box_11(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +Master_Triangle(p12,p32,s13,m12,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)
					 -Master_Triangle(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m32,mu2)
					 -Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_12(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +Master_Triangle(p12,p32,s13,m12,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)
					 -Master_Triangle(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m32,mu2)
					 -Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_13(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +Master_Triangle(p12,p32,s13,m12,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f2*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)
					 -Master_Triangle(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +Master_Triangle(p12,p22,s12,m12,m22,m32,mu2)
					 -Master_Triangle(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_211(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		    const double&  s12, const double&  s23,
		    const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       +PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2)
				       -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_212(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       +PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2)
				       -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_213(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       +PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2)
				       -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f2*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_11(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_11(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_221(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		    const double&  s12, const double&  s23,
		    const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_222(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_223(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_11(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f2*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_11(p12,p32,s13,m12,m32,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_231(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		    const double&  s12, const double&  s23,
		    const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_12(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_232(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_12(p22,p32,s23,m22,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_233(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		   const double&  s12, const double&  s23,
		   const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       +PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2)
				       -PV_Triangle_12(p22,p32,s23,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f2*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_12(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_12(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_12(p12,p22,s12,m12,m22,m42,mu2)
					 -2.*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_200(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		    const double&  s12, const double&  s23,
		    const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/(D-3.)*(2.*m12*Master_Box(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
		     +Master_Triangle(p22,p32,s23,m22,m32,m42,mu2)
		     -f1*PV_Box_11(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
		     -f2*PV_Box_12(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
		     -f3*PV_Box_13(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2));
}

METOOLS::DivArrC
METOOLS::PV_Box_3111(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
				       -4.*PV_Box_3001(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3112(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
				       -4.*PV_Box_3001(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3113(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -Master_Triangle(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
				       -4.*PV_Box_3001(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p12*p2p3)*(f2*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_211(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_21(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_21(p12,p22,s12,m12,m22,m42,mu2)));
}


METOOLS::DivArrC
METOOLS::PV_Box_3122(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_21(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
					 -4.*PV_Box_3002(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3123(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_223(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_23(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_23(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_223(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_23(p12,p32,s13,m12,m32,m42,mu2)
					 -2.*PV_Box_3003(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_223(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -2.*PV_Box_3002(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}


METOOLS::DivArrC
METOOLS::PV_Box_3222(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_21(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
					 -4.*PV_Box_3002(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3223(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_21(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p12*p2p3)*(f2*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_21(p12,p32,s13,m12,m32,m42,mu2)
					 -4.*PV_Box_3002(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_222(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3133(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_22(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -4.*PV_Box_3003(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3233(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_22(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -4.*PV_Box_3003(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3333(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_22(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p12*p2p3)*(f2*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_22(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_233(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 -4.*PV_Box_3003(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)));
}


METOOLS::DivArrC
METOOLS::PV_Box_3001(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p22*p32-sqr(p2p3))*(f1*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_24(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p3*p2p3-p1p2*p32)*(f2*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_24(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p2p3-p1p3*p22)*(f3*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3002(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p3*p2p3-p1p2*p32)*(f1*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_24(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p32-sqr(p1p3))*(f2*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_24(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p2p3*p12)*(f3*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)));
}

METOOLS::DivArrC
METOOLS::PV_Box_3003(const double&  p12, const double&  p22, const double&  p32, const double&  p42,
		     const double&  s12, const double&  s23,
		     const Complex& m12, const Complex& m22, const Complex& m32, const Complex& m42,
		     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  double p1p2 = 0.5*(s12-p12-p22);
  double p2p3 = 0.5*(s23-p22-p32);
  double s13 = p42-p22-2.*p1p2-2.*p2p3;
  double p1p3 = 0.5*(p42-p12-p22-p32-2.*p1p2-2.*p2p3);
  double det = p12*p22*p32-p12*sqr(p2p3)-p32*sqr(p1p2)-p22*sqr(p1p3)+2.*p1p2*p1p3*p2p3;
  Complex f1 = m22-m12-p12;
  Complex f2 = m32-m22-p22-2.*p1p2;
  Complex f3 = m42-m32-p32-2.*p2p3-2.*p1p3;
  return 0.5/det*((p1p2*p2p3-p1p3*p22)*(f1*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
				       -PV_Triangle_24(p22,p32,s23,m22,m32,m42,mu2)
				       +PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p1p2*p1p3-p12*p2p3)*(f2*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 +PV_Triangle_24(p12,p22,s12,m12,m22,m42,mu2)
					 -PV_Triangle_24(p12,p32,s13,m12,m32,m42,mu2))
		  +(p12*p22-sqr(p1p2))*(f3*PV_Box_200(p12,p22,p32,p42,s12,s23,m12,m22,m32,m42,mu2)
					 -PV_Triangle_22(p12,p22,s12,m12,m22,m42,mu2)
					 +PV_Triangle_22(p12,p22,s12,m12,m22,m32,mu2)));
}
