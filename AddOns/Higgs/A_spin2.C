
#include "ATOOLS/Math/MathTools.H"

using namespace ATOOLS;

/* =============================================================
   The resonant gg -> X -> gamma gamma  amplitudes, where X is spin-2 graviton,
   written in spinor products with phase given by specific definition
   of angle bracket spa(...) and square bracket spb(...)
   ============================================================= */

Complex ggXgamgam_mpmp(int i1, int i2, int i3, int i4, double kg=1.) {
  return -kg*pow(spa(i1,i3)*spb(i2,i4),2);
}

Complex ggXgamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3, int i4, int h4, double kg=1.) {
  if ( h1==h2 || h3==h4 ) return 0;
  else if ( h1==h3 ) {
    if ( h1==1 ) return ggXgamgam_mpmp(i2,i1,i4,i3,kg);
    else return ggXgamgam_mpmp(i1,i2,i3,i4,kg);
  }
  else { // h1==h4
    if ( h1==1 ) return ggXgamgam_mpmp(i2,i1,i3,i4,kg);
    else return ggXgamgam_mpmp(i1,i2,i4,i3,kg);
  }
}

//======

Complex qqXgamgam_mpmp(int i1, int i2, int i3, int i4, double kq=1.) {
  return -kq*pow(spa(i1,i3),2)*spb(i2,i4)*spb(i1,i4);
}
Complex qqXgamgam_pmpm(int i1, int i2, int i3, int i4, double kq=1.) {
  return -kq*pow(spb(i1,i3),2)*spa(i2,i4)*spa(i1,i4);
}

Complex qqXgamgam_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4, double kq=1.) {
  if ( h3==h4 ) return 0;
  else if ( h1==h3 ) {
    if ( h1==1 ) return qqXgamgam_pmpm(i1,i2,i3,i4,kq);
    else return qqXgamgam_mpmp(i1,i2,i3,i4,kq);
  }
  else { // h1==h4
    if ( h1==1 ) return qqXgamgam_pmpm(i1,i2,i4,i3,kq);
    else return qqXgamgam_mpmp(i1,i2,i4,i3,kq);
  }
}

//======

Complex qqgXgamgam_mppmp(int i1, int i2, int i3, int i4, int i5, double kg=1., double kq=1.) {
  return -kg*pow(spa(i1,i4),3)*pow(spb(i4,i5),2)*spa(i2,i4)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1)
         -(kq-kg)*pow(spa(i1,i4),2)/spa(i2,i3)/spa(i3,i1)*
          (spb(i5,i1)*spa(i1,i2)*spb(i2,i5)+spb(i5,i1)*spa(i1,i3)*spb(i3,i5)+spb(i5,i3)*spa(i3,i2)*spb(i2,i5));
}
Complex qqgXgamgam_pmmpm(int i1, int i2, int i3, int i4, int i5, double kg=1., double kq=1.) {
  return  kg*pow(spb(i1,i4),3)*pow(spa(i4,i5),2)*spb(i2,i4)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1)
         +(kq-kg)*pow(spb(i1,i4),2)/spb(i2,i3)/spb(i3,i1)*
          (spa(i5,i1)*spb(i1,i2)*spa(i2,i5)+spa(i5,i1)*spb(i1,i3)*spa(i3,i5)+spa(i5,i3)*spb(i3,i2)*spa(i2,i5));
}
Complex qqgXgamgam_pmpmp(int i1, int i2, int i3, int i4, int i5, double kg=1., double kq=1.) {
  return qqgXgamgam_mppmp(i2,i1,i3,i4,i5,kg,kq); }
Complex qqgXgamgam_mpmpm(int i1, int i2, int i3, int i4, int i5, double kg=1., double kq=1.) {
  return qqgXgamgam_pmmpm(i2,i1,i3,i4,i5,kg,kq); }

Complex qqgXgamgam_gen(int i1, int h1, int i2, int i3, int h3, int i4, int h4, int i5, int h5, double kg=1., double kq=1.) {
  if ( h4==h5 ) return 0;
  if ( h4==1 ) {
    if ( h1==1 ) {
      if ( h3==1 ) return qqgXgamgam_pmpmp(i1,i2,i3,i5,i4,kg,kq);
      else return qqgXgamgam_pmmpm(i1,i2,i3,i4,i5,kg,kq);
    }
    else { // h1==-1
      if ( h3==-1 ) return qqgXgamgam_mpmpm(i1,i2,i3,i4,i5,kg,kq);
      else return qqgXgamgam_mppmp(i1,i2,i3,i5,i4,kg,kq);
    }
  }
  else { // h4==-1
    if ( h1==1 ) {
      if ( h3==1 ) return qqgXgamgam_pmpmp(i1,i2,i3,i4,i5,kg,kq);
      else return qqgXgamgam_pmmpm(i1,i2,i3,i5,i4,kg,kq);
    }
    else { // h1==-1
      if ( h3==-1 ) return qqgXgamgam_mpmpm(i1,i2,i3,i5,i4,kg,kq);
      else return qqgXgamgam_mppmp(i1,i2,i3,i4,i5,kg,kq);
    }
  }
}

//======

Complex gggXgamgam_mppmp(int i1, int i2, int i3, int i4, int i5, double kg=1.) {
  return kg*pow(spa(i1,i4),4)*pow(spb(i4,i5),2)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1);
}
Complex gggXgamgam_pmmpm(int i1, int i2, int i3, int i4, int i5, double kg=1.) {
  return -kg*pow(spb(i1,i4),4)*pow(spa(i4,i5),2)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1);
}

Complex gggXgamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3, int i4, int h4, int i5, int h5, double kg=1.) {
  if ( h4==h5 ) return 0;
  if ( h1==h2 && h2==h3 ) return 0;
  else if ( h4==1 ) {
    if ( h1==1 ) {
      if ( h2==1 ) return gggXgamgam_mppmp(i3,i1,i2,i5,i4,kg);
      else if ( h3==1 ) return gggXgamgam_mppmp(i2,i3,i1,i5,i4,kg);
      else return gggXgamgam_pmmpm(i1,i2,i3,i4,i5,kg);
    }
    else { // h1==-1
      if ( h2==-1 ) return gggXgamgam_pmmpm(i3,i1,i2,i4,i5,kg);
      else if ( h3==-1 ) return gggXgamgam_pmmpm(i2,i3,i1,i4,i5,kg);
      else return gggXgamgam_mppmp(i1,i2,i3,i5,i4,kg);
    }
  }
  else { // h4==-1
    if ( h1==1 ) {
      if ( h2==1 ) return gggXgamgam_mppmp(i3,i1,i2,i4,i5,kg);
      else if ( h3==1 ) return gggXgamgam_mppmp(i2,i3,i1,i4,i5,kg);
      else return gggXgamgam_pmmpm(i1,i2,i3,i5,i4,kg);
    }
    else { // h1==-1
      if ( h2==-1 ) return gggXgamgam_pmmpm(i3,i1,i2,i5,i4,kg);
      else if ( h3==-1 ) return gggXgamgam_pmmpm(i2,i3,i1,i5,i4,kg);
      else return gggXgamgam_mppmp(i1,i2,i3,i4,i5,kg);
    }
  }
}

//=================================================================
// Special versions with fewer labels:

// g(1) g(2) -> gam(3) gam(4)
Complex ggXgamgam(int h1, int h2, int h3, int h4, double kg=1.) { return ggXgamgam_gen(1,h1,2,h2,3,h3,4,h4,kg); }

// q(1) qbar(2) -> gam(3) gam(4)
Complex qqbXgamgam(int h1, int h3, int h4, double kq=1.) { return qqXgamgam_gen(1,h1,2,3,h3,4,h4,kq); }
// qbar(1) q(2) -> gam(3) gam(4)
Complex qbqXgamgam(int h2, int h3, int h4, double kq=1.) { return qqXgamgam_gen(2,h2,1,3,h3,4,h4,kq); }

// g(1) g(2) -> gam(3) gam(4) g(5)
Complex ggXgamgamg(int h1, int h2, int h3, int h4, int h5, double kg=1.) { return gggXgamgam_gen(1,h1,2,h2,5,h5,3,h3,4,h4,kg); }

// q(1) qbar(2) -> gam(3) gam(4) g(5)
Complex qqbXgamgamg(int h1, int h3, int h4, int h5, double kg=1., double kq=1.) { return qqgXgamgam_gen(1,h1,2,5,h5,3,h3,4,h4,kg,kq); }
// qbar(1) q(2) -> gam(3) gam(4) g(5)
Complex qbqXgamgamg(int h2, int h3, int h4, int h5, double kg=1., double kq=1.) { return qqgXgamgam_gen(2,h2,1,5,h5,3,h3,4,h4,kg,kq); }
// q(1) g(2) -> gam(3) gam(4) q(5)
Complex qgXgamgamq(int h1, int h2, int h3, int h4, double kg=1., double kq=1.) { return qqgXgamgam_gen(1,h1,5,2,h2,3,h3,4,h4,kg,kq); }
// g(1) q(2) -> gam(3) gam(4) q(5)
Complex gqXgamgamq(int h1, int h2, int h3, int h4, double kg=1., double kq=1.) { return qqgXgamgam_gen(2,h2,5,1,h1,3,h3,4,h4,kg,kq); }

//========

// running couplings for kg and kq
double kgr(double muR, double kg=1., double kq=1., double muR0=0.) {
  if (kq==kg) return kg;
  double lmu2=(muR0<=0.)?2.*log(Flavour(kf_h0).Mass()/muR):2.*log(muR0/muR);
  double Z1g=alpha_s(muR)/8./PI*8.0/3.0*N_f/2.0;
  double Z1q=alpha_s(muR)/8./PI*16.0/3.0*4.0/3.0;
  return (exp((Z1g+Z1q)*lmu2)*(kg-kq)*Z1g+kq*Z1g+kg*Z1q)/(Z1g+Z1q);
}

double kqr(double muR, double kg=1., double kq=1., double muR0=0.) {
  if (kq==kg) return kq;
  double lmu2=(muR0<=0.)?2.*log(Flavour(kf_h0).Mass()/muR):2.*log(muR0/muR);
  double Z1g=alpha_s(muR)/8./PI*8.0/3.0*N_f/2.0;
  double Z1q=alpha_s(muR)/8./PI*16.0/3.0*4.0/3.0;
  return (exp((Z1g+Z1q)*lmu2)*(kq-kg)*Z1q+kq*Z1g+kg*Z1q)/(Z1g+Z1q);
}

// g(1) g(2) -> gam(3) gam(4) 1-loop
Complex ggXgamgam1l(int h1, int h2, int h3, int h4, double muR, double kg=1., double kq=1.) {
  if (kq==kg) return 0.;
  Complex L12 = log(muR*muR/sij(1,2))+(sij(1,2)>0.0?I*PI:0.);
  return ggXgamgam(h1,h2,h3,h4,kg) * 1.0/4.0*L12*(8.0/3.0*N_f/2.0)*(kg-kq); }

// q(1) qbar(2) -> gam(3) gam(4) 1-loop
Complex qqbXgamgam1l(int h1, int h3, int h4, double muR, double kg=1., double kq=1.) {
  if (kq==kg) return 0.;
  double mh = Flavour(kf_h0).Mass();
  Complex L12 = log(muR*muR/sij(1,2))+(sij(1,2)>0.0?I*PI:0.);
  return qqbXgamgam(h1,h3,h4,kq) * 1.0/4.0*L12*(16.0/3.0*4.0/3.0)*(kq-kg); }
// qbar(1) q(2) -> gam(3) gam(4) 1-loop
Complex qbqXgamgam1l(int h2, int h3, int h4, double muR, double kg=1., double kq=1.) {
  if (kq==kg) return 0.;
  Complex L12 = log(muR*muR/sij(1,2))+(sij(1,2)>0.0?I*PI:0.);
  return qbqXgamgam(h2,h3,h4,kq) * 1.0/4.0*L12*(16.0/3.0*4.0/3.0)*(kq-kg); }

