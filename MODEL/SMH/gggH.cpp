#include "gggH.h"
#include "gpl.h"
#include <cmath>
#include <iostream>

const double z2 = M_PI*M_PI/6.;

using namespace std;

void Fggg(double mb, double mh, double s, double t, double u,
	  std::complex<double> *F) {
  double mh2 = mh*mh;
  std::complex<double> T1=Omegagggppp1l(mb, mh, s, t, u);//+++
  std::complex<double> T2=Omegagggppm1l(mb, mh, s, t, u);//++-
  std::complex<double> T3=Omegagggppm1l(mb, mh, t, u, s);//-++
  std::complex<double> T4=Omegagggppm1l(mb, mh, u, s, t);//+-+
  F[0]+=mb*mb*(mh2*mh2*T1+s*s*T2-t*t*T3-u*u*T4)/(2.*mh2*s*u);
  F[1]+=mb*mb*(mh2*mh2*T1-s*s*T2+t*t*T3-u*u*T4)/(2.*mh2*s*t);
  F[2]+=mb*mb*(mh2*mh2*T1-s*s*T2-t*t*T3+u*u*T4)/(2.*mh2*t*u);
  F[3]+=mb*mb*(-mh2*mh2*T1+s*s*T2+t*t*T3+u*u*T4)/(mh2*s*t*u);
}

complex<double> F1ggg(double mb, double mh, double s12, double s13, double s23) {
  double mh2 = mh*mh;
  return 0.5*(mh2/s12/s23*Omegagggppp1l(mb, mh, s12, s13, s23) + 1./mh2*(s12/s23*Omegagggppm1l(mb, mh, s12, s13, s23) + s13*s13/s12/s23*Omegagggppm1l(mb, mh, s13, s12, s23) + s23/s12*Omegagggppm1l(mb, mh, s23, s13, s12)));
}

complex<double> F2ggg(double mb, double mh, double s12, double s13, double s23) {
  return -s23/s12*F1ggg(mb, mh, s13, s12, s23) - s23/2.*(F4ggg(mb, mh, s13, s12, s23) - F4ggg(mb, mh, s12, s13, s23));
}

complex<double> F3ggg(double mb, double mh, double s12, double s13, double s23) {
  return -s12/s13*F1ggg(mb, mh, s23, s13, s12) - s12/2.*(F4ggg(mb, mh, s23, s13, s12) - F4ggg(mb, mh, s12, s13, s23));
}

// TODO: Check the relative sign in front of the last term...
complex<double> F4ggg(double mb, double mh, double s12, double s13, double s23) {
  double mh2 = mh*mh;
  return -mh2/s12/s13/s23*Omegagggppp1l(mb, mh, s12, s13, s23) - 1./mh2*(s12/s13/s23*Omegagggppm1l(mb, mh, s12, s13, s23) + s13/s12/s23*Omegagggppm1l(mb, mh, s13, s12, s23) - s23/s12/s13*Omegagggppm1l(mb, mh, s23, s13, s12));
}

complex<double> Omegagggppp1l(double mb, double mh, double s12, double s13, double s23) {
  if (s12 > 0 && s13 > 0 && s23 > 0 && s13 < mh*mh && s23 < mh*mh && s12 <= mh*mh - s13 - s23) {
    return Omegagggppp1l1a(mb, mh, s12, s13, s23);
  } else if (s12 > 0 && s13 < 0 && s23 < 0) {
    return Omegagggppp1l2a(mb, mh, s12, s13, s23);
  } else if (s12 < 0 && s13 > 0 && s23 < 0) {
    // There is no region 3a here, since we deal with 3 identical particles, just swap the mandelstams s13 and s23
    return -Omegagggppp1l(mb, mh, s12, s23, s13);
  } else if (s12 < 0 && s13 < 0 && s23 > 0) {
    return Omegagggppp1l4a(mb, mh, s12, s13, s23);
  } else {
    return 0.0;
  }
}

complex<double> Omegagggppm1l(double mb, double mh, double s12, double s13, double s23) {
  if (s12 > 0 && s13 > 0 && s23 > 0 && s13 < mh*mh && s23 < mh*mh && s12 <= mh*mh - s13 - s23) {
    return Omegagggppm1l1a(mb, mh, s12, s13, s23);
  } else if (s12 > 0 && s13 < 0 && s23 < 0) {
    return Omegagggppm1l2a(mb, mh, s12, s13, s23);
  } else if (s12 < 0 && s13 > 0 && s23 < 0) {
    // There is no region 3a here, since we deal with 3 identical particles, just swap the mandelstams s13 and s23
    return -Omegagggppm1l(mb, mh, s12, s23, s13);
  } else if (s12 < 0 && s13 < 0 && s23 > 0) {
    return Omegagggppm1l4a(mb, mh, s12, s13, s23);
  } else {
    return 0.0;
  }
}



complex<double> Omegagggppp1l1a(double mb, double mh, double s12, double s13, double s23) {
  double y = s13/(mh*mh);
  double z = s23/(mh*mh);
  double lnbH = log(mb*mb/(mh*mh));
  double lnz = log(z);
  double lny = log(y);
  double ln1mz = log(1. - z);
  double ln1yz = log(1. + y/z);
  double ln1my1mz = log(1. - y/(1. - z));
  
  double re = -4. + ln1my1mz*ln1my1mz/2. + ln1mz*ln1mz/2. + (3.*lnbH*lnbH)/2. + lny*lny/2. + lnbH*(-lny - lnz) + ln1mz*(-2.*ln1yz - lnbH + lny - lnz) + lny*lnz + lnz*lnz/2.
                  + ln1my1mz*(ln1mz - lnbH + lny + lnz) - 8.*z2 - 2.*GPL(1.,0.,y) - 2.*GPL(-z,1. - z,y);
  double im = -ln1my1mz - ln1mz + 3.*lnbH - lny - lnz;
  
  return complex<double>(re,M_PI*im);
}

complex<double> Omegagggppp1l2a(double mb, double mh, double s12, double s13, double s23) {
  double u1 = -s13/s12;
  double v1 = (mh*mh)/s12;
  double lnbH = log(mb*mb/(mh*mh));
  double lnu = log(u1);
  double lnv = log(v1);
  double ln1mu = log(1. - u1);
  double ln1mv = log(1. - v1);
  double ln1uv = log(1. + u1/v1);
  double ln1mu1mv = log(1. - u1/(1. - v1));
  
  double re = -4. + ln1mu1mv*ln1mu1mv/2. + ln1mv*ln1mv/2. + (3.*lnbH*lnbH)/2. + lnu*lnu/2. + ln1mu1mv*(ln1mv - lnbH + lnu - 3.*lnv) + ln1mv*(-2.*ln1mu - lnbH + lnu - lnv)
                  + 2.*ln1mu*lnv + 2.*ln1uv*lnv - 3.*lnu*lnv + (5.*lnv*lnv)/2. + lnbH*(-lnu + 3.*lnv) - 8.*z2 - 2.*GPL(1.,1. - v1,u1) - 2.*GPL(-v1,0.,u1);
  double im = -2.*ln1mu + ln1mu1mv + ln1mv - 2.*ln1uv + lnbH + lnu - lnv;
  
  return complex<double>(re,M_PI*im);
}

complex<double> Omegagggppp1l4a(double mb, double mh, double s12, double s13, double s23) {
  double u2 = -s13/s23;
  double v2 = (mh*mh)/s23;
  double lnbH = log(mb*mb/(mh*mh));
  double lnu = log(u2);
  double lnv = log(v2);
  double ln1mu = log(1. - u2);
  double ln1mv = log(1. - v2);
  double ln1uv = log(1. + u2/v2);
  double ln1mu1mv = log(1. - u2/(1. - v2));
  
  double re = -4. + ln1mu1mv*ln1mu1mv/2. + ln1mv*ln1mv/2. + (3.*lnbH*lnbH)/2. + lnu*lnu/2. + ln1mu1mv*(ln1mv - lnbH + lnu - 3.*lnv) + ln1mv*(-2.*ln1mu - lnbH + lnu - lnv)
                  + 2.*ln1mu*lnv + 2.*ln1uv*lnv - 3.*lnu*lnv + (5.*lnv*lnv)/2. + lnbH*(-lnu + 3.*lnv) - 8.*z2 - 2.*GPL(1.,1. - v2,u2) - 2.*GPL(-v2,0.,u2);
  double im = -2.*ln1mu + ln1mu1mv + ln1mv - 2.*ln1uv + lnbH + lnu - lnv;
  
  return complex<double>(re,M_PI*im);
}


complex<double> Omegagggppm1l1a(double mb, double mh, double s12, double s13, double s23) {
  double y = s13/(mh*mh);
  double z = s23/(mh*mh);
  
  double lnbH = log(mb*mb/(mh*mh));
  double lnz = log(z);
  double lny = log(y);
  double ln1mz = log(1. - z);
  double ln1yz = log(1. + y/z);
  double ln1my1mz = log(1. - y/(1. - z));
  
  double denyzm1 = 1./(y + z - 1.);
  double denz = 1./(z - 1.);
  double deny = 1./(y - 1.);
  
  double re = denyzm1*denyzm1*(-4.*lny*y*z*deny*deny*(-2. + 2.*y + z) - 4.*lnz*y*z*denz*denz*(-2. + y + 2.*z))
            + denyzm1*denyzm1*denyzm1*(4.*y*z*GPL(1.,0.,y) + (lny*lnz - GPL(1.,0.,z))*(1. - 2.*y + y*y - 2.*z - 2.*y*z + z*z) + 4.*z2*(1. - 2.*y + y*y - 2.*z + y*z + z*z))
            + denyzm1*(-ln1my1mz*ln1my1mz/2. - ln1my1mz*ln1mz - ln1mz*ln1mz/2. + 2.*ln1mz*ln1yz + ln1my1mz*lnbH + ln1mz*lnbH - lnbH*lnbH/2. - ln1my1mz*lny - ln1mz*lny - ln1my1mz*lnz
                               + GPL(0.,1.,z) + 2.*GPL(-z,1. - z,y) + (-(lnbH*lnz) + lnz*lnz/2.)*denz*(-1. + 2.*y + z) + deny*((-(lnbH*lny) + lny*lny/2.)*(-1. + y + 2.*z) + 4.*denz*(1. - 2.*y + y*y - 2.*z + y*z + z*z)));
  double im = denyzm1*(ln1my1mz + ln1mz - lnbH - lnz*denz*(-1. + 2.*y + z) - lny*deny*(-1. + 2.*z + y));
  
  return complex<double>(re,M_PI*im);
}



complex<double> Omegagggppm1l2a(double mb, double mh, double s12, double s13, double s23) {
  double u1 = -s13/s12;
  double v1 = (mh*mh)/s12;
  double lnbH = log(mb*mb/(mh*mh));
  double denu = 1./(u1 - 1.);
  double denv = 1./(v1 - 1.);
  double denuv = 1./(u1 + v1);
  double lnu = log(u1);
  double lnv = log(v1);
  double ln1mu = log(1. - u1);
  double ln1mv = log(1. - v1);
  double ln1uv = log(1. + u1/v1);
  double ln1mu1mv = log(1. - u1/(1. - v1));

  double re = -4.*lnu*u1*v1*denuv*denuv*(-1. + u1 + v1)*(1. + u1 + v1) + v1*(lnbH*lnbH/2. + u1*(4.*ln1mu*ln1mv - 4.*ln1mu*lnv - 4.*ln1uv*lnv + 4.*GPL(1., 1. - v1, u1) + 4.*GPL(-v1, 0., u1))*(-1. + u1 + v1)
            - ln1mu1mv*lnu*(1. - 4.*u1 + 4.*u1*u1 + 4.*u1*v1) - ln1mv*lnu*(1. - 4.*u1 + 4.*u1*u1 + 4.*u1*v1)) + denu*denu*((4.*ln1mu1mv + 4.*ln1mv)*u1*v1*(-2. + u1)*(-1. + u1 + v1)
            - 4.*lnv*u1*v1*denuv*denuv*(-1. + u1 + v1)*(-1. + u1 - u1*u1 - v1 - 2.*u1*v1 + u1*u1*v1 - 2.*v1*v1 + u1*v1*v1)) + v1*denuv*((-(lnbH*lnu) + lnu*lnu/2.)*(-2. + u1 + v1)
            + lnu*lnv*(2. - u1 - 4.*u1*u1 + 4.*u1*u1*u1 - v1 - 4.*u1*v1 + 8.*u1*u1*v1 + 4.*u1*v1*v1)) + denu*(v1*((ln1mu1mv*ln1mu1mv/2. + ln1mu1mv*ln1mv + ln1mv*ln1mv/2. - ln1mu1mv*lnbH - ln1mv*lnbH - GPL(0., 1., v1))*(1. + u1)
            + ln1mu1mv*lnv*(-1. + 3.*u1 - 8.*u1*u1 + 4.*u1*u1*u1 - 4.*u1*v1 + 4.*u1*u1*v1) + GPL(1., 0., v1)*(-3. + 5.*u1 - 8.*u1*u1 + 4.*u1*u1*u1 - 4.*u1*v1 + 4.*u1*u1*v1)) + v1*denuv*(4.*(1. - u1 + u1*u1 + u1*v1)
            + lnbH*lnv*(2. - 3.*u1 + 3.*u1*u1 - v1 + 3.*u1*v1) - (lnv*lnv*(-2. + 3.*u1 + u1*u1 - 8.*u1*u1*u1 + 4.*u1*u1*u1*u1 + v1 + u1*v1 - 12.*u1*u1*v1 + 8.*u1*u1*u1*v1 - 4.*u1*v1*v1 + 4.*u1*u1*v1*v1))/2.
            + 2.*z2*(3. - 3.*u1 + 11.*u1*u1 - 16.*u1*u1*u1 + 8.*u1*u1*u1*u1 + 11.*u1*v1 - 24.*u1*u1*v1 + 16.*u1*u1*u1*v1 - 8.*u1*v1*v1 + 8.*u1*u1*v1*v1)));
  double im = -(lnbH*v1*denu*denuv*(2. - u1 + u1*u1 + v1 + u1*v1)) + (-ln1mu1mv - ln1mv - lnu + lnv)*v1*(1. - 4.*u1 + 4.*u1*u1 + 4.*u1*v1)
              + u1*v1*(-1. + u1 + v1)*(4.*ln1mu + 4.*ln1uv + 4.*denu*denu*denuv*denuv*(-1. + u1 - u1*u1 - v1 - 2.*u1*v1 + u1*u1*v1 - 2.*v1*v1 + u1*v1*v1));
  
  return complex<double>(re,M_PI*im);     
}


complex<double> Omegagggppm1l4a(double mb, double mh, double s12, double s13, double s23) {
  double u2 = -s13/s23;
  double v2 = (mh*mh)/s23;
  double lnbH = log(mb*mb/(mh*mh));
  double denv = 1./(v2 - 1.);
  double denuv = 1./(u2 + v2);
  double denuv1 = 1./(u2 + v2 - 1.);
  double lnu = log(u2);
  double lnv = log(v2);
  double ln1mu = log(1. - u2);
  double ln1mv = log(1. - v2);
  double ln1uv = log(1. + u2/v2);
  double ln1mu1mv = log(1. - u2/(1. - v2));
  
  double re = v2*(denuv1*(ln1mu1mv*ln1mu1mv/2. - 2.*ln1mu*ln1mv + ln1mu1mv*ln1mv + ln1mv*ln1mv/2. - ln1mu1mv*lnbH - ln1mv*lnbH + lnbH*lnbH/2. + ln1mu1mv*lnu + ln1mv*lnu + 2.*ln1mu*lnv - 3.*ln1mu1mv*lnv
            - GPL(0., 1., v2) - 2.*GPL(1., 1. - v2, u2) + denuv*((lnbH*lnu - lnu*lnu/2.)*(-2. + u2 + v2) + denv*(-4.*(1. - u2 + u2*u2 - 2.*v2 + 2.*u2*v2 + v2*v2) - lnbH*lnv*(2. - u2 + 2.*u2*u2 - 3.*v2 + 3.*u2*v2 + v2*v2))))
            + denuv1*denuv1*denuv*denuv*(-4.*lnu*u2*(-1. + 2.*u2 + 2.*v2) + 4.*lnv*u2*denv*denv*(-1. + 2.*u2 - 2.*u2*u2 + u2*u2*u2 + 4.*v2 - 8.*u2*v2 + 4.*u2*u2*v2 - 7.*v2*v2 + 7.*u2*v2*v2 + 4.*v2*v2*v2))
            + denuv1*denuv1*denuv1*(u2*(-4.*ln1uv*lnv + 4.*GPL(-v2, 0., u2)) - GPL(1., 0., v2)*(3. - 2.*u2 + 3.*u2*u2 - 6.*v2 + 6.*u2*v2 + 3.*v2*v2)
            + denuv*(lnu*lnv*(-2. + 5.*u2 + u2*u2*u2 + 5.*v2 - 4.*u2*v2 + 3.*u2*u2*v2 - 4.*v2*v2 + 3.*u2*v2*v2 + v2*v2*v2)
            - 2.*z2*(-3. + 9.*u2 - 11.*u2*u2 + 3.*u2*u2*u2 + 9.*v2 - 20.*u2*v2 + 9.*u2*u2*v2 - 9.*v2*v2 + 9.*u2*v2*v2 + 3.*v2*v2*v2)
            - (lnv*lnv*denv*(2. - 5.*u2 + 2.*u2*u2 - 5.*u2*u2*u2 + 2.*u2*u2*u2*u2 - 7.*v2 + 11.*u2*v2 - 11.*u2*u2*v2 + 7.*u2*u2*u2*v2 + 9.*v2*v2 - 11.*u2*v2*v2 + 9.*u2*u2*v2*v2 - 5.*v2*v2*v2 + 5.*u2*v2*v2*v2 + v2*v2*v2*v2))/2.)));
  double im = v2*(denuv1*(-2.*ln1mu + ln1mu1mv + ln1mv + lnu + lnbH*denuv*(-2. + u2 + v2)) - 4.*u2*denuv*denuv*denuv1*denuv1*(-1. + 2.*u2 + 2.*v2)
            + denuv1*denuv1*denuv1*(4.*ln1uv*u2 - lnv*denv*(-1. + 8.*u2 - 5.*u2*u2 + 2.*u2*u2*u2 + 3.*v2 - 12.*u2*v2 + 5.*u2*u2*v2 - 3.*v2*v2 + 4.*u2*v2*v2 + v2*v2*v2)));
  
  return complex<double>(re,M_PI*im);     
}
