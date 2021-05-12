#include "qqgH.h"
#include "gpl.h"
#include <cmath>
#include <iostream>

const double z2 = M_PI*M_PI/6.;

using namespace std;

complex<double> F1qqg(double mb, double mh, double s13, double s23) {
  double s12 = mh*mh - s13 - s23;
  return s12*mh*mh*Omegaqqgmpp1l(mb, mh, s12, s13, s23);
}

complex<double> F2qqg(double mb, double mh, double s13, double s23) {
  double s12 = mh*mh - s13 - s23;
  return s12*mh*mh*Omegaqqgmpp1l(mb, mh, s12, s23, s13);
}

complex<double> Omegaqqgmpp1l(double mb, double mh, double s12, double s13, double s23) {
  if (s12 > 0 && s13 > 0 && s23 > 0 && s13 < mh*mh && s23 < mh*mh && s12 <= mh*mh - s13 - s23) {
    cout << "1a" << endl;
    return Omegaqqgmpp1l1a(mb, mh, s12, s13, s23);
  } else if (s12 > 0 && s13 < 0 && s23 < 0) {
    cout << "2a" << endl;
    return Omegaqqgmpp1l2a(mb, mh, s12, s13, s23);
  } else if (s12 < 0 && s13 > 0 && s23 < 0) {
    cout << "3a" << endl;
    return Omegaqqgmpp1l3a(mb, mh, s12, s13, s23);
  } else if (s12 < 0 && s13 < 0 && s23 > 0) {
    cout << "4a" << endl;
    return Omegaqqgmpp1l4a(mb, mh, s12, s13, s23);
  } else {
    return 0.0;
  }
}

complex<double> Omegaqqgmpp1l1a(double mb, double mh, double s12, double s13, double s23) {
  double y = s13/(mh*mh);
  double z = s23/(mh*mh);
  double lnbH = log(mb*mb/(mh*mh));
  double denyz = 1./(y + z);
  double lnz = log(z);
  double lny = log(y);
  double ln1mz = log(1. - z);
  double ln1yz = log(1. + y/z);
  double ln1my1mz = log(1. - y/(1. - z));
  
  double re = (-4. - ln1my1mz*ln1my1mz - 2.*ln1my1mz*ln1mz - ln1mz*ln1mz + 2.*ln1my1mz*lnbH + 2.*ln1mz*lnbH)*denyz + (4.*ln1my1mz + 4.*ln1mz)*denyz*denyz*(-1. + y + z);
  double im = (2.*ln1my1mz + 2.*ln1mz)*denyz;
  
  return complex<double>(re,M_PI*im);
}

complex<double> Omegaqqgmpp1l2a(double mb, double mh, double s12, double s13, double s23) {
  double v2 = (mh*mh)/s12;
  double lnbH = log(mb*mb/(mh*mh));
  double denv = 1./(v2 - 1.);
  double lnv = log(v2);
  
  double re = v2*((-4. - 2.*lnbH*lnv - lnv*lnv)*denv + 4.*lnv*denv*denv);
  double im = -2.*lnv*v2*denv;
  
  return complex<double>(re,M_PI*im);
}

complex<double> Omegaqqgmpp1l3a(double mb, double mh, double s12, double s13, double s23) {
  double u3 = -s23/s13;
  double v3 = (mh*mh)/s13;
  double lnbH = log(mb*mb/(mh*mh));
  double denu = 1./(u3 - 1.);
  double lnv = log(v3);
  double ln1mv = log(1. - v3);
  double ln1mu1mv = log(1. - u3/(1. - v3));
  
  double re = v3*((4. + ln1mu1mv*ln1mu1mv + 2.*ln1mu1mv*ln1mv + ln1mv*ln1mv - 2.*ln1mu1mv*lnbH - 2.*ln1mv*lnbH - 2.*ln1mu1mv*lnv - 2.*ln1mv*lnv + 2.*lnbH*lnv + lnv*lnv + 6.*z2)*denu
            + (-4.*ln1mu1mv - 4.*ln1mv + 4.*lnv)*denu*denu*(-1. + u3 + v3));
  double im = -2.*lnbH*v3*denu - 4.*v3*denu*denu*(-1. + u3 + v3);
  
  return complex<double>(re,M_PI*im);
}

complex<double> Omegaqqgmpp1l4a(double mb, double mh, double s12, double s13, double s23) {
  double u4 = -s13/s23;
  double v4 = (mh*mh)/s23;
  double lnbH = log(mb*mb/(mh*mh));
  double denu = 1./(u4 - 1.);
  double lnv = log(v4);
  double ln1mv = log(1. - v4);
  double ln1mu1mv = log(1. - u4/(1. - v4));
  
  double re = v4*((4. + ln1mu1mv*ln1mu1mv + 2.*ln1mu1mv*ln1mv + ln1mv*ln1mv - 2.*ln1mu1mv*lnbH - 2.*ln1mv*lnbH - 2.*ln1mu1mv*lnv - 2.*ln1mv*lnv + 2.*lnbH*lnv + lnv*lnv + 6.*z2)*denu
            + (-4.*ln1mu1mv - 4.*ln1mv + 4.*lnv)*denu*denu*(-1. + u4 + v4));
  double im = -2.*lnbH*v4*denu - 4.*v4*denu*denu*(-1. + u4 + v4);
  
  return complex<double>(re,M_PI*im);
}
