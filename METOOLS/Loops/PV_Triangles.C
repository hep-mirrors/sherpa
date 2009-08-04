#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! C_10(p12,p22,s12;m02,m12,m22)
//! C_11(p12,p22,s12;m02,m12,m22)
//! C_12(p12,p22,s12;m02,m12,m22)
//! C_20(p2,p2,0;m02,m2,m2)
//! C_21(p12,p22,s12;m02,m12,m22)
//! C_22(p12,p22,s12;m02,m12,m22)
//! C_23(p12,p22,s12;m02,m12,m22)
//! C_24(p2,p2,0;m02,m2,m2)
//! C_24(p12,p22,s12;m02,m12,m22)
//! C_30(p2,p2,0;m02,m2,m2)
//! C_38(p2,p2,0;m02,m2,m2)
/*! all momenta incoming

        p2
       -----|\              p1 + p2 + p3 = 0
            | \ m0
            |  \            s12 = p3^2 = (-p1-p2)^2 = (p1+p2)^2
            |   \   p3
         m2 |    |-----
            |   /
            |  /            s12 > 0    s-channel
            | / m1          s12 < 0    t-channel
       -----|/              s12 = 0    in wave-function corrections
        p1
*/
METOOLS::DivArrC
METOOLS::PV_Triangle_10(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (!IsZero(s12) || !IsEqual(p12,p22) || !IsEqual(m12,m22)) {
    THROW(fatal_error,"call in ill-defined situation");
  }
  else {
    double p2(0.5*(p12+p22));
    Complex m2(0.5*(m12+m22));
    if (IsZero(p2)) {
      //! C_10(0,0,0;m02,m2,m2) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else {
      //! C_10(p2,p2,0;m02,m2,m2) = 1/(2*p2) * ( B_0(p2;m12,m2) - B_0(0;m2,m2)
      //!                                         + fC C_0(p2,p2,0;m02,m2,m2) )
      Complex fC(-p2+m2-m02);
      return 0.5/p2*(Master_Bubble(p2,m02,m2,mu2)
                     -Master_Bubble(0.,m2,m2,mu2)
                     +(!IsZero(fC)?fC*Master_Triangle(p2,p2,0.,m02,m2,m2,mu2)
                                  :DivArrC(0.,0.,0.,0.,0.,0.)));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_11(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (IsZero(s12)) {
    //! s12 = 0  ==>  p1 = p2
    //!   ==> C_11(p12,p22,0;m02,m12,m22) = C_10(p2,p2,0;m02,m12,m22)
    return PV_Triangle_10(p12,p22,s12,m02,m12,m22,mu2);
  }
  else {
    //! C_11(p12,p22,s12;m02,m12,m22)
    //!   = 1/2 * 1/(p12p22-(p1p2)^2)
    //!        * ( p22 ( B_0(p22;m02,m22) - B_0(s12;m12,m22)
    //!                   + f1C C_0(p12,p22,s12;m02,m12,m22) )
    //!           - (p1p2) ( B_0(p12;m02,m12) - B_0(s;m12,m22)
    //!                       + f2C C_0(p12,p22,s12;m02,m12,m22) ) )
    double p1p2(0.5*(s12-p12-p22));
    Complex f1C(-p12+m12-m02), f2C(-p22+m22-m02);
    if (IsZero(p12*p22-sqr(p1p2))) {
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else {
      DivArrC C0(((!IsZero(f1C)) || (!IsZero(f2C)))?
                      Master_Triangle(p12,p22,s12,m02,m12,m22,mu2):
                      DivArrC(0.,0.,0.,0.,0.,0.));
      return 0.5/(p12*p22-sqr(p1p2))
             * ((!IsZero(p22)?p22*(Master_Bubble(p22,m02,m22,mu2)
                                   -Master_Bubble(s12,m12,m22,mu2)
                                   +f1C*C0)
                             :DivArrC(0.,0.,0.,0.,0.,0.))
                -(!IsZero(p1p2)?p1p2*(Master_Bubble(p12,m02,m12,mu2)
                                      -Master_Bubble(s12,m12,m22,mu2)
                                      +f2C*C0)
                               :DivArrC(0.,0.,0.,0.,0.,0.)));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_12(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  //! C_12(p12,p22,s12;m02,m12,m22) = C_11(p22,p12,s12;m02,m22,m12)
  return PV_Triangle_11(p22,p12,s12,m02,m22,m12,mu2);
}


METOOLS::DivArrC
METOOLS::PV_Triangle_20(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (!IsZero(s12) || !IsEqual(p12,p22) || !IsEqual(m12,m22)) {
    THROW(fatal_error,"call in ill-defined situation");
  }
  else {
    double p2(0.5*(p12+p22));
    if (IsZero(p2)) {
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else {
      Complex m2(0.5*(m12+m22));
      Complex fC(-p2+m2-m02);
      //! C_20(p2,p2,0;m02,m2,m2)
      //!   = 1/2 * 1/(p4(D-1))
      //!        * ( D/2 (p2 B_1(p2;m02,m2) + p2 B_0(0;m2,m2)
      //!                   + fC p2 C_10(p2,p2,0;m02,m2,m2) )
      //!           - p2 ( B_0(0;m2,m2)
      //!                  + m02 C_0(p2,p2,0;m02,m2,m2)))
      return 0.5/(p2*(D-1.))
              *(0.5*D*(PV_Bubble_1(p2,m02,m2,mu2)
                       +Master_Bubble(0.,m2,m2,mu2)
                       +(!IsZero(fC)?fC*PV_Triangle_10(p2,p2,0.,m02,m2,m2,mu2)
                                    :DivArrC(0.,0.,0.,0.,0.,0.)))
                - (Master_Bubble(0,m2,m2,mu2)
                   +(!IsZero(m02)?m02*Master_Triangle(p2,p2,0,m02,m2,m2,mu2)
                                 :DivArrC(0.,0.,0.,0.,0.,0.))));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_21(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! C_21(p12,p22,s12;m02,m12,m22)
  //!   = 1/2 * 1/(p12p22-(p1p2)^2)
  //!        * ( p22 ( B_1(s12;m12,m22) - B_0(s12;m12,m22)
  //!                   + f1C C_11(p12,p22,s12;m02,m12,m22)
  //!                   - 2 C_24(p12,p22,s12;m02,m12,m22) )
  //!           - (p1p2) ( B_1(p12;m02,m12) - B_1(s12;m22,m12)
  //!                       + f2C C_12(p12,p22,s12;m02,m12,m22) ) )
  double p1p2(0.5*(s12-p12-p22));
  if (IsZero(p12*p22-sqr(p1p2))) {
    msg_Out()<<"not implemented yet\n";
    return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  else {
    Complex f1C(-p12+m12-m02), f2C(-p22+m22-m02);
    DivArrC C11(((!IsZero(f1C)) || (!IsZero(f2C)))?
                    PV_Triangle_24(p12,p22,s12,m02,m12,m22,mu2):
                    DivArrC(0.,0.,0.,0.,0.,0.));
    return 0.5/(p12*p22-sqr(p1p2))
            *((!IsZero(p22)?p22*(PV_Bubble_1(s12,m12,m22,mu2)
                                 -Master_Bubble(s12,m12,m22,mu2)
                                 +f1C*C11
                                 -2.*PV_Triangle_24(p12,p22,s12,m02,m12,m22,mu2))
                           :DivArrC(0.,0.,0.,0.,0.,0.))
              -(!IsZero(p1p2)?p1p2*(PV_Bubble_1(p12,m02,m12,mu2)
                                    -PV_Bubble_1(s12,m22,m12,mu2)
                                    +f2C*C11)
                             :DivArrC(0.,0.,0.,0.,0.,0.)));
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_22(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! C_22(p12,p22,s12;m02,m12,m22) = C_21(p22,p12,s12;m02,m22,m12)
  return PV_Triangle_21(p22,p12,s12,m02,m22,m12,mu2);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_23(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! C_23(p12,p22,s12;m02,m12,m22)
  //!   = 1/2 * 1/(p12p22-(p1p2)^2)
  //!        * ( p12 ( B_1(p12;m02,m12) - B_1(s12;m22,m12)
  //!                   + f2C C_12(p12,p22,s12;m02,m12,m22) )
  //!           - (p1p2) ( B_1(s12;m12,m12) + B_0(s12;m12,m22)
  //!                       + f1C C_11(p12,p22,s12;m02,m12,m22)
  //!                      - 2 C_24(p12,p22,s12;m02,m12,m22) ) )
  double p1p2(0.5*(s12-p12-p22));
  if (IsZero(p12*p22-sqr(p1p2))) {
    msg_Out()<<"not implemented yet\n";
    return DivArrC(0.,0.,0.,0.,0.,0.);
  }
  else {
    Complex f1C(-p12+m12-m02), f2C(-p22+m22-m02);
    DivArrC C11(((!IsZero(f1C)) || (!IsZero(f2C)))?
                    PV_Triangle_24(p12,p22,s12,m02,m12,m22,mu2):
                    DivArrC(0.,0.,0.,0.,0.,0.));
    return 0.5/(p12*p22-sqr(p1p2))
            *((!IsZero(p12)?p12*(PV_Bubble_1(p12,m02,m12,mu2)
                                 -PV_Bubble_1(s12,m22,m12,mu2)
                                 +f2C*C11)
                           :DivArrC(0.,0.,0.,0.,0.,0.))
              -(!IsZero(p1p2)?p1p2*(PV_Bubble_1(s12,m12,m22,mu2)
                                    -Master_Bubble(s12,m12,m22,mu2)
                                    +f1C*C11
                                    -2.*PV_Triangle_24(p12,p22,s12,m02,m12,m22,mu2))
                             :DivArrC(0.,0.,0.,0.,0.,0.)));        
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_24(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (IsZero(s12) && IsEqual(p12,p22) && IsEqual(m12,m22)) {
    //! C_24(p2,p2,0;m02,m2,m2)
    //!   = 1/(p4(D-1)) * ( p4 ( B_0(0;m2,m2) + m02 C_0(p2,p2,0;m02,m2,m2) )
    //!                    - p2/2 ( p2 B_1(p2;m02,m2) + p2 B_0(0;m2,m2)
    //!                             + fC p2 C_10(p2,p2,0;m02,m2,m2) ) )
    double p2(0.5*(p12+p22));
    Complex m2(0.5*(m12+m22));
    Complex fC(-p2+m2-m02);
    return 1./(D-1.)*(Master_Bubble(0.,m2,m2,mu2)
                      +(!IsZero(m02)?m02*Master_Triangle(p2,p2,0.,m02,m2,m2,mu2)
                                    :DivArrC(0.,0.,0.,0.,0.,0.))
                      -0.5*(PV_Bubble_1(p2,m02,m2,mu2)
                            +Master_Bubble(0.,m2,m2,mu2)
                            +(!IsZero(fC)?fC*PV_Triangle_10(p2,p2,0.,m02,m2,m2,mu2)
                                         :DivArrC(0.,0.,0.,0.,0.,0.))));
  }
  else {
    //! C_24(p12,p22,s12;m02,m12,m22)
    //!   = 1/(D-2) * ( m02 C_0(p12,p22,s12;m02,m12,m22)
    //!                 + 1/2 B_0(s12;m12,m22)
    //!                 - f1C C_11(p12,p22,s12;m02,m12,m22)
    //!                 - f2C C_12(p12,p22,s12;m02,m12,m22) )
    Complex f1C(-p12+m12-m02), f2C(-p22+m22-m02);
    return 1./(D-2.)*((!IsZero(m02)?m02*Master_Triangle(p12,p22,s12,m02,m12,m22,mu2)
                                   :DivArrC(0.,0.,0.,0.,0.,0.))
                      +0.5*Master_Bubble(s12,m12,m22,mu2)
                      -(!IsZero(f1C)?f1C*PV_Triangle_11(p12,p22,s12,m02,m12,m22,mu2)
                                    :DivArrC(0.,0.,0.,0.,0.,0.))
                      -(!IsZero(f2C)?f2C*PV_Triangle_12(p12,p22,s12,m02,m12,m22,mu2)
                                    :DivArrC(0.,0.,0.,0.,0.,0.)));
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}


METOOLS::DivArrC
METOOLS::PV_Triangle_30(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (!IsZero(s12) || !IsEqual(p12,p22) || !IsEqual(m12,m22)) {
    THROW(fatal_error,"call in ill-defined situation");
  }
  else {
    double p2(0.5*(p12+p22));
    if (IsZero(p2)) {
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else {
      //! C_30(p2,p2,0;m02,m2,m2)
      //!   = 1/(2p2) * ( B_21(p2;m02,m2) - B_21(0;m2,m2)
      //!                   + fC C_20(p2,p2,0;m02,m2,m2) )
      //!     - 2/p2 C_38(p2,p2,0;m02,m2,m2)
      //!
      //!   = 1/(2p2) * ( B_21(p2;m02,m2) - B_21(0;m2,m2)
      //!                   + fC C_20(p2,p2,0;m02,m2,m2) )
      //!     - 1/p4 * ( B_22(p2;m02,m2) - B_22(0;m2,m2)
      //!                  + fC C_24(p2,p2,0;m02,m2,m2) )
      Complex m2(0.5*(m12+m22));
      Complex fC(-p2+m2-m02);
      return 0.5/p2*(PV_Bubble_21(p2,m02,m2,mu2)
                     -PV_Bubble_21(0.,m2,m2,mu2)
                     +(!IsZero(fC)?fC*PV_Triangle_20(p2,p2,0,m02,m2,m2,mu2)
                                  :DivArrC(0.,0.,0.,0.,0.,0.)))
             -1./sqr(p2)*(PV_Bubble_22(p2,m02,m2,mu2)
                          -PV_Bubble_22(0.,m2,m2,mu2)
                          +(!IsZero(fC)?fC*PV_Triangle_24(p2,p2,0,m02,m2,m2,mu2)
                                  :DivArrC(0.,0.,0.,0.,0.,0.)));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}

METOOLS::DivArrC
METOOLS::PV_Triangle_38(const double&  p12, const double&  p22, const double&  s12,
                        const Complex& m02, const Complex& m12, const Complex& m22,
                        double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (!IsZero(s12) || !IsEqual(p12,p22) || !IsEqual(m12,m22)) {
    THROW(fatal_error,"call in ill-defined situation");
  }
  else {
    double p2(0.5*(p12+p22));
    if (IsZero(p2)) {
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else {
      //! C_38(p2,p2,0;m02,m2,m2)
      //!   = 1/(2p2) * ( B_22(p2;m02,m2) - B_22(0;m2,m2)
      //!                   + fC C_24(p2,p2,0;m02,m2,m2) )
      Complex m2(0.5*(m12+m22));
      Complex fC(-p2+m2-m02);
      return 0.5/sqr(p2)*(PV_Bubble_22(p2,m02,m2,mu2)
                          -PV_Bubble_22(0.,m2,m2,mu2)
                          +(!IsZero(fC)?fC*PV_Triangle_24(p2,p2,0,m02,m2,m2,mu2)
                                  :DivArrC(0.,0.,0.,0.,0.,0.)));

    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}




