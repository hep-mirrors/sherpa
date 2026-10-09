/*!
  \file ee_aa_QED_Poles.C

  See ee_aa_QED_Poles.H. Reference: yfs-helpers/diphoton-two-loop,
  massive_poles.py (mpmath), compared by cmp_poles.py to <= 1e-12.
*/
#include "EXTRA_XS/NNLO/ee_aa_QED_Poles.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H"
#include <cmath>

namespace {

  typedef std::complex<double> C;
  const double s_pi  = 3.14159265358979323846264338327950288;
  const double s_pi2 = s_pi*s_pi;
  const double s_z2  = s_pi2/6.;
  const double s_z3  = 1.20205690315959428539973816151144999;

  //! c = a*b for series a (lowest power la, na terms), b (lb, nb), keeping
  //! powers la+lb .. kmax
  template <class T>
  void Mul(const T *a, int la, int na, const T *b, int lb, int nb, T *c, int kmax)
  {
    const int lc(la + lb);
    for (int k(lc); k <= kmax; ++k) c[k - lc] = T(0.);
    for (int i(0); i < na; ++i)
      for (int j(0); j < nb; ++j) {
        const int k(la + i + lb + j);
        if (k <= kmax) c[k - lc] += a[i]*b[j];
      }
  }

  //! Catani's I^(1) for e+e- gamma gamma, mu^2 = s, N_F' = 0, eps^-2 .. eps^2 (AGT eq. 3.6)
  void I1Series(C i1[5])
  {
    i1[0] = C(-1., 0.);
    i1[1] = C(-1.5, -s_pi);
    i1[2] = C(7.*s_pi2/12., -1.5*s_pi);
    i1[3] = C(s_z3/3. + 7.*s_pi2/8., s_pi*s_pi2/4.);
    i1[4] = C(-121.*s_pi2*s_pi2/1440. + 0.5*s_z3, s_pi*s_z3/3. + 3.*s_pi*s_pi2/8.);
  }

  //! Mitov-Moch Z^(1) (eps^-2 .. eps^2) and the C_F^2 part of Z^(2) (eps^-4 .. eps^0), alpha/4pi
  void ZSeries(double L, double z1[5], double z2[5])
  {
    const double L2(L*L), L3(L2*L), L4(L3*L);
    z1[0] = 2.; z1[1] = 2.*L + 1.; z1[2] = L2 + L + 4. + s_z2;
    z1[3] = L3/3. + L2/2. + (4. + s_z2)*L + 8. + s_z2/2. - 2.*s_z3/3.;
    z1[4] = L4/12. + L3/6. + (2. + s_z2/2.)*L2 + (8. + s_z2/2. - 2.*s_z3/3.)*L
      + 16. + 2.*s_z2 - s_z3/3. + 9./20.*s_z2*s_z2;
    z2[0] = 2.; z2[1] = 4.*L + 2.; z2[2] = 4.*L2 + 4.*L + 17./2. + 2.*s_z2;
    z2[3] = 8./3.*L3 + 4.*L2 + (17. + 4.*s_z2)*L + 83./4. - 4.*s_z2 + 32./3.*s_z3;
    z2[4] = 4./3.*L4 + 8./3.*L3 + (17. + 4.*s_z2)*L2 + (83./2. - 8.*s_z2 + 64./3.*s_z3)*L
      + 561./8. + 61./2.*s_z2 - 22./3.*s_z3 - 48.*std::log(2.)*s_z2 - 77./5.*s_z2*s_z2;
  }

}

namespace EXTRAXS { namespace ee_aa_QED_Poles {

void TOverT0(double x, double t[3])
{
  const double R(x/(1. - x) + (1. - x)/x), Q(1./(x*(1. - x)));
  t[0] = 1.; t[1] = -1. - Q/R; t[2] = Q/R;
}

void OneLoopSeries(double x, std::complex<double> m[5])
{
  using ee_aa_QED_TwoLoop::Li;
  const double y(1. - x), X(std::log(x)), Y(std::log(y));
  const double Li2x(Li(2, x)), Li2y(Li(2, y)), Li3x(Li(3, x)), Li3y(Li(3, y)),
    Li4x(Li(4, x)), Li4y(Li(4, y)), Li4z(Li(4, (x - 1.)/x)), Li4zp(Li(4, -x/y));
  const double z3(s_z3);
  (void)Li4x;
#include "ee_aa_QED_MSeries.inc"
}

void MasslessSeries(double x, double n02[5], double n11[5], double m1[5])
{
  C m[5], i1[5], i1c[5], tmp[9];
  OneLoopSeries(x, m);
  I1Series(i1);
  for (int k(0); k < 5; ++k) { m1[k] = 2.*m[k].real(); i1c[k] = std::conj(i1[k]); }
  // I2 = -I1^2/2 + e^{eps g}/Gamma(1-eps) H/(4 eps), through eps^0 (eps^-4 .. eps^0)
  const double H(2.*(s_pi2/2. - 6.*s_z3 - 3./8.));
  C i2[5];
  Mul(i1, -2, 5, i1, -2, 5, tmp, 0);
  for (int k(0); k < 5; ++k) i2[k] = -0.5*tmp[k];
  i2[3] += C(H/4., 0.);              // eps^-1; the eps^0 term of the normalisation is 0
  // P02/T = 2 Re(I1 m + I2), P11/T = -|I1|^2 + 2 Re(I1^* m), eps^-4 .. eps^0
  C im[5], ic[5], ii[5];
  Mul(i1, -2, 5, m, -2, 5, im, 0);
  Mul(i1c, -2, 5, m, -2, 5, ic, 0);
  Mul(i1, -2, 5, i1c, -2, 5, ii, 0);
  for (int k(0); k < 5; ++k) {
    n02[k] = 2.*(im[k] + i2[k]).real();
    n11[k] = -ii[k].real() + 2.*ic[k].real();
  }
  n02[4] += ee_aa_QED_TwoLoop::Fin02(x);
  n11[4] += ee_aa_QED_TwoLoop::Fin11(x);
}

void MassiveSeries(double x, double L, double m1m[3], double m2m[3], bool ol)
{
  double n02[5], n11[5], m1[5], z1[5], z2[5], t[9];
  MasslessSeries(x, n02, n11, m1);
  ZSeries(L, z1, z2);
  // m1m = m1 + Z1 (Z1 = 2 z1, alpha/2pi), eps^-2 .. eps^1 (eps^-2 cancels)
  for (int k(0); k < 3; ++k) m1m[k] = m1[k + 1] + z1[k + 1];
  // m2m = m2 + Z1 m1 + Z1^2/4 + Z2/2, eps^-4 .. eps^0 (eps^-4, eps^-3 cancel)
  double a[5], b[5];
  Mul(z1, -2, 5, m1, -2, 5, a, 0);
  Mul(z1, -2, 5, z1, -2, 5, b, 0);
  for (int k(0); k < 5; ++k) t[k] = n02[k] + n11[k] + a[k] + 0.25*b[k] + 0.5*z2[k];
  for (int k(0); k < 3; ++k) m2m[k] = t[k + 2];
  if (ol) {
    // R = 1 + zeta2/2 eps^2 per loop: shifts only these orders
    m1m[2] += 0.5*s_z2*m1m[0];
    m2m[2] += s_z2*m2m[0];
  }
}


} }
