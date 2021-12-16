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
  if (IsZero(s) || s < pow(10.,-6.)) {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_1(0;0,0) = ?
      return -0.5*Master_Bubble(s,0.,0.,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_1(0;0,m2) = -0.5*B_0(0;0,m2) + 1/4
      return -0.5*Master_Bubble(s,m02,m12,mu2) + 0.25*DivArrC(0.,0.,0.,1.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      //! B_1(0;m2,m2) = -0.5*B_0(0;m2,m2)
      return -0.5*Master_Bubble(s,m02,m12,mu2);
    }
    else if (IsEqual(s,m02) && abs(m12/m02) > pow(10.,6.)) {
      Complex fb(-s+m12-m02);
      return 0.5/s*(Master_Tadpole(m02,mu2)
		    -Master_Tadpole(m12,mu2)
		    +fb*Master_Bubble(s,m02,m12,mu2));
    }
    else {
      //! B_1(0;m12,m22) = -0.5*B_0(0;m12,m22) + 0.5*(m12-m22)*B_0p(0;m12,m22)
      return -0.5*Master_Bubble(s,m02,m12,mu2) + 0.5*(m02-m12)*Master_Bubble_Prime(s,m02,m12,mu2);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_1(s;0,0) = ?
      return -0.5*Master_Bubble(s,0.,0.,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsZero(m02)) {
	if (IsEqual(s,m2)) {
	  //! B_1(m2;0,m2) = -1/2*A_0(m2)/m2
	  return -0.5*Master_Tadpole(m2,mu2)/m2;
	}
	else {
	  //! B_1(s;0,m2) = ?
	  return 0.5/s*(-Master_Tadpole(m2,mu2)+
			(m2-s)*Master_Bubble(s,0.,m2,mu2));
	}
      }
      else if (IsZero(m12)) { // Note: B_1(p2,m02,m12) = - B_1(p2,m12,m02) - B_0(p2,m02,m12)
	if (IsEqual(s,m2)) {
	  //! B_1(m2;m2,0) = 1/2*A_0(m2)/m2
	  return 0.5*Master_Tadpole(m2,mu2)/m2-Master_Bubble(s,0.,m2,mu2);
	}
	else {
	  //! B_1(s;m2,0) = ?
	  return 0.5/s*(Master_Tadpole(m2,mu2)-
			(m2+s)*Master_Bubble(s,0.,m2,mu2));
	}
      }
    }	
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_1(m2;m2,m2) = ?
        msg_Out()<<"B_1 all equal not implemented yet\n";
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

// Derivative of B_1
METOOLS::DivArrC
METOOLS::PV_Bubble_1_Prime(const double& s,
			   const Complex& m02, const Complex& m12,
			   double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  //! two massless internal lines
  if (IsZero(s)) {
    //! B_1_prime(0;m02,m12) = ?
    if (IsZero(m12) && IsZero(m02)) {
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      return DivArrC(0.,0.,0.,-1./(12.*m12),0.,0.);
    }
    else if (IsZero(m12)) {
      return DivArrC(0.,0.,0.,-1./(3.*m02),0.,0.);
    }
    else if (IsZero(m02)) {
      return DivArrC(0.,0.,0.,-1./(6.*m12),0.,0.);
    }  
  }
  if (IsZero(m02) && IsZero(m12)) {
    //! B_1_prime(s;0,0) = 1/(2 s)
    return DivArrC(0.,0.,0.,0.5/s,0.,0.);
  }
  if (IsZero(m12) && IsEqual(s,m02)) {
    //! B_1_prime(s;m2,0) = 1/(2 m2)*(1/epsIR + 3. + ln(mu2/m02)
    return 0.5/m02*DivArrC(0.,1.,0.,3.+log(mu2/m02),0.,0.);
  }
  if (IsZero(m02) && IsEqual(s,m12)) {
    return DivArrC(0.,0.,0.,-0.5/s,0.,0.);
  }
  if (IsZero(m02) && !IsEqual(s,m12)) {
    // Direct calculation Bardin/Passarino eq. 5.41
    if (abs(s/m12) < pow(10.,-8.)) { return DivArrC(0.,0.,0.,-1./(6.*m12)-s/(12.*sqr(m12)),0.,0.); }
    return -1./(2.*pow(s,3.))*DivArrC(0.,0.,0.,s*(s-2.*m12)+2.*m12*(s-m12)*CLog(1.-s/m12,-1),0.,0.);
  }
  if (IsEqual(m02,m12)) {
    return -0.5*Master_Bubble_Prime(s,m02,m12,mu2);
  }
  if (IsEqual(s,m02) && !IsEqual(m02,m12)) {
    Complex xm, xp;
    if (abs(m02/m12) > pow(10.,-6.)) {
      xm = 1./(2.*m02)*(-sqrt(4.*sqr(m02)+sqr(m12))+m12);
    }
    else {
      xm = -m02/m12;
    }
    xp = 1./(2.*m02)*(sqrt(4.*sqr(m02)+sqr(m12))+m12);
    // Bardin/Passarino eq. 5.43 and series 5.45
    if (abs(xp) > pow(10.,9.)) {
      return DivArrC(0.,0.,0.,1./(m02*(xp-xm))*(1./(3.*xp)-1./(3.*xm)-1./(4.*xp)+1./(4.*xm)),0.,0.);
    }
    // Bardin/Passarino eq. 5.43
    return DivArrC(0.,0.,0.,1./(m02*(xp-xm))
      *(-sqr(xp)*CLog(1.-1./xp)-xp-1./2.+sqr(xm)*CLog(1.-1./xm)+xm+1./2.
  	+pow(xp,3.)*CLog(1.-1./xp)+sqr(xp)+xp/2.+1./3.-pow(xm,3.)*CLog(1.-1./xm)-sqr(xm)-xm/2.-1./3.),
  		   0.,0.);
  }
  else {
    //! B_1_prime(s,m02,m12) = -(m12-m02)/(2.*s^2)*(B_0(s,m02,m12)-B_0(0,m02,m12)) + (m12-m02-s)/(2*s)*B_0_prime(s,m02,m12)
    if (IsZero(m02) && abs(m12/s) > pow(10.,6.)) return -0.25/s-1./(12.*m12)+0.5*(m12-s)/s*Master_Bubble_Prime(s,m02,m12,mu2);
    return -0.5*(m12-m02)/sqr(s)*(Master_Bubble(s,m02,m12,mu2)-Master_Bubble(0.,m02,m12,mu2)) + 0.5*(m12-m02-s)/s*Master_Bubble_Prime(s,m02,m12,mu2);
  }
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
      return 1./3.*Master_Bubble(s,m02,m12,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_21(0;0,m2) = 1./3.*B_0(0;0,m2) - 4/9
      return 1./3.*Master_Bubble(s,m02,m12,mu2)-4./9.*DivArrC(0.,0.,0.,1.,0.,0.);
    }
    else if (IsEqual(m02,m12)) {
      //! B_21(0;m2,m2) = 1./3.*B_0(0;m2,m2)
      return 1./3.*Master_Bubble(s,m02,m12,mu2);
    }
    else {
      //! B_21(0;m12,m22) = ?
      msg_Out()<<"B_21 2 masses not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_21(s;0,0) = -D/(2*(D-1))*B_1(s,0,0)
      return -0.5*D/(D-1)*PV_Bubble_1(s,m02,m12,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_21(m2;0,m2) = 1/2*(D-2)/(D-1)*A_0(m2)/m2
        return 0.5*(D-2.)/(D-1.)*Master_Tadpole(m2,mu2)/m2;
      }
      else {
        //! B_21(s;0,m2) = 1/(s*(D-1)) * ( 1/2*(D-2)*A_0(m2)
        //!                              + D/2*(m2-s)*B_1(s;0,m2))
        return 1./(s*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
			      + 0.5*D*(m12 - s)*PV_Bubble_1(s,m02,m12,mu2));
      }
    }
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_21(m2;m2,m2) = ?
        msg_Out()<<"B_21 m2 m2 m2 not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_21(s;m2,m2) = 1/(s*(D-1)) * ( 1/2*(D-2)*A_0(m2)
        //!                                - D/2*s*B_1(s;m2,m2)
	//!                                - m2*B_0(s,m2,m2) )
        return 1./(s*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
			      - 0.5*D*s*PV_Bubble_1(s,m02,m12,mu2)
			      - m02*Master_Bubble(s,m02,m12,mu2));
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
METOOLS::PV_Bubble_21_Prime(const double& s,
			    const Complex& m02, const Complex& m12,
			    double mu2=0.) {
  if (mu2 == 0.)   mu2 = GLOBAL_RENORMALISATION_SCALE;
  if (IsZero(s)) {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_21(0;0,0) = ?
      msg_Out()<<"B_21 all 0 not implemented yet\n";
      return 1./3.*Master_Bubble_Prime(s,m02,m12,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      //! B_21(0;0,m2) = 1./3.*B_0p(0;0,m2)
      msg_Out()<<"B_21 one mass not implemented yet\n";
      return 1./3.*Master_Bubble_Prime(s,m02,m12,mu2);
    }
    else if (IsEqual(m02,m12)) {
      //! B_21(0;m2,m2) = 1./3.*B_0p(0;m2,m2)
      return 1./3.*Master_Bubble_Prime(s,m02,m12,mu2);
    }
    else {
      //! B_21(0;m12,m22) = ?
      msg_Out()<<"B_21 2 masses not implemented yet\n";
      return DivArrC(0.,0.,0.,0.,0.,0.);
    }
  }
  else {
    if (IsZero(m02) && IsZero(m12)) {
      //! B_21(s;0,0) = -D/(2*(D-1))*B_1(s,0,0)
      return -0.5*D/(D-1)*PV_Bubble_1_Prime(s,m02,m12,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_21(m2;0,m2) = 1/2*(D-2)/(D-1)*A_0(m2)/m2
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_21(s;0,m2) = -1/(s^2*(D-1)) * ( 1/2*(D-2)*A_0(m2)
        //!                                   + D/2*m2*B_1(s;0,m2))
	//!                -1/(s*(D-1)) * D/2*B_1p
        return -1./(sqr(s)*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
				    + 0.5*D*m12*PV_Bubble_1(s,m02,m12,mu2))
	  -1./(sqr(s)*(D-1.))*0.5*D*PV_Bubble_1_Prime(s,m02,m12,mu2);
      }
    }
    else if (IsEqual(m02,m12)) {
      if (IsEqual(s,0.5*(m02+m12))) {
        //! B_21(m2;m2,m2) = ?
        msg_Out()<<"B_21 m2 m2 m2 not implemented yet\n";
        return DivArrC(0.,0.,0.,0.,0.,0.);
      }
      else {
        //! B_21(s;m2,m2) = -1/(s^2*(D-1)) * ( 1/2*(D-2)*A_0(m2)
	//!                                   - m2*B_0(s,m2,m2) )
        //!                 +1/(s*(D-1)) * (- D/2*s*B_1p(s;m2,m2)
	//!                                 - m2*B_0p(s,m2,m2) )
        return -1./(sqr(s)*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
				    - m02*Master_Bubble(s,m02,m12,mu2))
	  +1./(s*(D-1.))*(- 0.5*D*s*PV_Bubble_1_Prime(s,m02,m12,mu2)
			  - m02*Master_Bubble_Prime(s,m02,m12,mu2));
      }
    }
    else {
      //! B_21(s;m12,m22) = 1/(s*(D-1)) * ( 1/2*(D-2)*A_0(m22)
      //!                                  + D/2*fb*B_1(s;m12,m22)
      //!                                  - m12*B_0(s;m12,m22) )
      Complex fb(-s+m12-m02);
      return -1./(sqr(s)*(D-1.))*(0.5*(D-2.)*Master_Tadpole(m12,mu2)
				  +0.5*D*(m12-m02)*PV_Bubble_1(s,m02,m12,mu2)
				  -m02*Master_Bubble(s,m02,m12,mu2))
	+1./(s*(D-1.))*(- 0.5*D*s*PV_Bubble_1_Prime(s,m02,m12,mu2)
			- m02*Master_Bubble_Prime(s,m02,m12,mu2));
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
     return -0.5*Master_Tadpole((m02+m12)/2.,mu2);
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
	return 0.5/(D-1.)*s*PV_Bubble_1(s,m02,m12,mu2);
    }
    else if (IsZero(m02) || IsZero(m12)) {
      Complex m2(m02+m12);
      if (IsEqual(s,m2)) {
        //! B_22(m2;0,m2) = 1/2*1/(D-1)*A_0(m2)
        return 0.5/(D-1.)*Master_Tadpole(m2,mu2);
      }
      else {
        //! B_22(s;0,m2) = ?
 	return 1./(D-1.)*(0.5*Master_Tadpole(m12,mu2)
			  +0.5*(s-m12)*PV_Bubble_1(s,m02,m12,mu2));
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
	return 1./(D-1.)*(0.5*Master_Tadpole(m12,mu2)
			  +0.5*s*PV_Bubble_1(s,m02,m12,mu2)
			  +m02*Master_Bubble(s,m02,m12,mu2));
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





