#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace std;
using namespace ATOOLS;
using namespace METOOLS;

//! B_1(s;m02,m12), B_21(s;m02,m12), B_22(s;m02,m12)
//! in Passarino-Veltman reduction
/*!
            -------
      q   /   m1   \   q
    -----|         |------    s = q^2
         \   m2   /
          -------
*/

METOOLS::DivArrC
METOOLS::PV_Bubble_1(const double& s,
                     const Complex& m02, const Complex& m12,
                     double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! two massless internal lines
  if (IsZero(s)) {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_1(0;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_1(0;0,m2) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      //! B_1(0;m2,m2) = -B_0(0;m2,m2)
      return -Master_Bubble(s,m02,m12,mu2);
    }
    else {
      //! B_1(0;m12,m22) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_1(s;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_1(m2;0,m2) = -1/2*A_0(m2)/m2
        return -0.5*Master_Tadpole(m2,mu2)/m2;
      }
      else {
        //! B_1(s;0,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
    }
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_1(m2;m2,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_1(s;m2,m2) = -1/2 * B_0(s,;m2,m2)
        return -0.5*Master_Bubble(s,m02,m12,mu2);
      }
    }
    else {
      //! B_1(s;m12,m22) = 1/2s * (A_0(m02)-A_0(m12)+fb*B_0(s;m02,m12)
      Complex fb(-s+m12-m02);
      return 0.5/s*(Master_Tadpole(m02,mu2)
                    -Master_Tadpole(m12,mu2)
                    +fb*Master_Bubble(s,m02,m12,mu2));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}


METOOLS::DivArrC
METOOLS::PV_Bubble_21(const double& s,
                      const Complex& m02, const Complex& m12,
                      double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! two massless internal lines
  if (IsZero(s)) {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_21(0;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_21(0;0,m2) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      //! B_21(0;m2,m2) = B_0(0;m2,m2)
      return -Master_Bubble(s,m02,m12,mu2);
    }
    else {
      //! B_21(0;m12,m22) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_21(s;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_21(m2;0,m2) = 1/2*(D-2)/(D-1)*A_0(m2)/m2
        return 0.5*(D-2.)/(D-1.)*Master_Tadpole(m2,mu2)/m2;
      }
      else {
        //! B_21(s;0,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
    }
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_21(m2;m2,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_21(s;m2,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
    }
    else {
      //! B_21(s;m12,m22) = 1/(s*(D-1)) * ( 1/2*(D-2)*A_0(m22)
      //!                                  + D/2*fb*B_1(s;m12,m22)
      //!                                  - m12*B_0(s;m12,m22) )
      Complex fb(-s+m12-m02);
      return 1./(s*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
                            +0.5*D*fb*PV_Bubble_1(s,m02,m12,mu2)
                            -m02*Master_Bubble(s,m02,m12,mu2));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}


METOOLS::DivArrC
METOOLS::PV_Bubble_22(const double& s,
                      const Complex& m02, const Complex& m12,
                      double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! two massless internal lines
  if (IsZero(s)) {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_22(0;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_22(0;0,m2) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      //! B_22(0;m2,m2) = 1/2*A_0(m2)
     return 0.5*Master_Tadpole((m02+m12)/2.,mu2);
    }
    else {
      //! B_22(0;m12,m22) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_22(s;0,0) = ?
      msg_Out()<<"not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_22(m2;0,m2) = 1/2*1/(D-1)*A_0(m2)
        return 0.5/(D-1.)*Master_Tadpole(m2,mu2);
      }
      else {
        //! B_22(s;0,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
    }
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_22(m2;m2,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_22(s;m2,m2) = ?
        msg_Out()<<"not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
    }
    else {
      //! B_22(s;m12,m22) = 1/(D-1) * ( 1/2*A_0(m22)
      //!                              - 1/2*fb*B_1(s;m12,m22)
      //!                              + m12*B_0(s;m12,m22) )
      Complex fb(-s+m12-m02);
      return 1./(D-1.)*(0.5*Master_Tadpole(m12,mu2)
                        -0.5*fb*PV_Bubble_1(s,m02,m12,mu2)
                        +m02*Master_Bubble(s,m02,m12,mu2));
    }
  }
  return DivArrC(0.,0.,0.,0.,0.,0.);
}





