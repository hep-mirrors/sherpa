/*!
  \file ee_aa_QED_TwoLoop.C

  Two-loop QED for e+e- -> gamma gamma (Anastasiou, Glover, Tejeda-Yeomans,
  hep-ph/0201274), see ee_aa_QED_TwoLoop.H. Transcription of App. A.1 (A_s,
  B_s, E_4;s, F_2;s) and App. B.1 (G_1;s, X_4;s, X_5;s); G_1;s and r1 are
  verified against Recola per-helicity amplitudes (1e-11, 1e-8.
*/
#include "EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H"
#include <cmath>
#include <stdexcept>

namespace {

  const double s_pi  = 3.14159265358979323846264338327950288;
  const double s_pi2 = s_pi*s_pi;
  const double s_pi4 = s_pi2*s_pi2;
  const double s_z2  = s_pi2/6.;
  const double s_z3  = 1.20205690315959428539973816151144999;
  const double s_z4  = s_pi4/90.;

  double Zeta(int k)
  {
    // zeta(k) for k = 2, 3, 4 and for k <= 0 (the Li_n expansion at w -> 1)
    switch (k) {
    case 4: return s_z4;
    case 3: return s_z3;
    case 2: return s_z2;
    case 0: return -0.5;
    default: break;
    }
    if (k > 0) throw std::runtime_error("Zeta: unexpected argument");
    // zeta(1 - 2m) = -B_2m/(2m), zeta(-2m) = 0
    static const double B2m[] = {1./6., -1./30., 1./42., -1./30., 5./66., -691./2730.,
                                 7./6., -3617./510., 43867./798., -174611./330.,
                                 854513./138., -236364091./2730., 8553103./6.,
                                 -23749461029./870., 8615841276005./14322.};
    if (k % 2 == 0) return 0.;              // zeta(-2m) = 0
    const int m((1 - k)/2);                  // k = 1 - 2m
    if (m - 1 >= int(sizeof(B2m)/sizeof(double))) return 0.;
    return -B2m[m - 1]/(2.*m);
  }

  double LiSeries(int n, double w)
  {
    double sum(0.), p(w);
    for (int k(1); k < 200; ++k) {
      const double term(p/std::pow(double(k), n));
      sum += term;
      if (std::abs(term) < 1e-18*std::abs(sum)) break;
      p *= w;
    }
    return sum;
  }

  //! Li_n(e^mu) for mu < 0, |mu| < 2 pi: sum_k zeta(n-k) mu^k/k! + log term
  double LiNearOne(int n, double w)
  {
    const double mu(std::log(w));
    double H(0.);
    for (int j(1); j < n; ++j) H += 1./j;
    double sum(0.), p(1.), fact(1.);
    // |mu| < ln 2 here: the terms fall like (mu/2 pi)^k; zeta(-2m) = 0 makes
    // every second term vanish, so no early exit on a zero term
    for (int k(0); k < n + 28; ++k) {
      if (k > 0) { p *= mu; fact *= k; }
      if (k == n - 1) sum += p/fact*(H - std::log(-mu));
      else sum += Zeta(n - k)*p/fact;
    }
    return sum;
  }

}

namespace EXTRAXS { namespace ee_aa_QED_TwoLoop {

double Li(int n, double w)
{
  if (n < 2 || n > 4) throw std::runtime_error("Li: n must be 2, 3 or 4");
  if (w > 1.) throw std::runtime_error("Li: real argument > 1");
  if (w == 0.) return 0.;
  if (w == 1.) return Zeta(n);
  if (std::abs(w) <= 0.5) return LiSeries(n, w);
  if (w > 0.5) return LiNearOne(n, w);
  if (w >= -1.) return std::pow(2., 1 - n)*Li(n, w*w) - Li(n, -w);
  // w < -1: inversion, L = ln(-w)
  const double L(std::log(-w)), iw(1./w);
  switch (n) {
  case 2: return -Li(2, iw) - s_pi2/6. - 0.5*L*L;
  case 3: return  Li(3, iw) - s_pi2/6.*L - L*L*L/6.;
  default: return -Li(4, iw) - 7.*s_pi4/360. - s_pi2/12.*L*L - L*L*L*L/24.;
  }
}

namespace {

  struct Kin {
    double x, y, X, Y, L2x, L2y, L3x, L3y, L4x, L4y, L4z, tu, ts2, tu2, tu3;
    explicit Kin(double xx)
    {
      x = xx; y = 1. - xx; X = std::log(x); Y = std::log(y);
      L2x = Li(2, x); L2y = Li(2, y); L3x = Li(3, x); L3y = Li(3, y);
      L4x = Li(4, x); L4y = Li(4, y); L4z = Li(4, (x - 1.)/x);
      tu = x/y; ts2 = x*x; tu2 = tu*tu; tu3 = tu2*tu;
    }
  };

  template <class G> double Sym(G g, double x, double S)
  { return g(Kin(x), S) + g(Kin(1. - x), S); }

  double gA(const Kin &k, double)
  {
    const double X(k.X), Y(k.Y);
    const double a(128.*k.L4z - 128.*k.L4x + 128.*k.L4y + (-64./3. + 128.*Y)*k.L3x
                   + (64./3.*X - 64./3.*s_pi2)*k.L2x + 16./3.*X*X*X*X - 64./3.*X*X*X*Y
                   + (-16. + 32./3.*s_pi2 + 32.*Y*Y)*X*X + (-64./3.*s_pi2*Y + 48. + 160./9.*s_pi2)*X
                   + 64./3.*s_z3 + 224./45.*s_pi4 - 128.*Y*s_z3);
    const double b(32./3.*k.L3x - 32./3.*k.L3y + (-32./3.*X - 32./3.*Y)*k.L2x
                   + (-32./3.*Y*Y - 80./9.*s_pi2 - 64./3.)*X + (64./3. + 32./3.*s_pi2)*Y);
    const double c(416./3.*k.L3x + 64.*k.L3y*X - 416./3.*k.L2x*X + (8.*Y*Y + 16.)*X*X
                   + (-8./3.*Y + 80./3. + 112./9.*s_pi2 - 64.*s_z3 - 64.*Y*Y)*X
                   - 416./3.*s_z3 - 148./9.*s_pi2 + 44./45.*s_pi4);
    return a*k.tu + b*k.ts2 + 24.*X*X*k.tu2 + c;
  }

  double gB(const Kin &k, double S)
  {
    const double X(k.X), Y(k.Y);
    const double a(-112.*k.L4z - 88.*k.L4y + (-128.*Y + 48.*X - 64.)*k.L3x
                   + (-16.*Y - 16.*X + 12.)*k.L3y
                   + (12.*Y - 4.*Y*Y + 8.*X*X - 8.*s_pi2 + 64.*X)*k.L2x
                   + 2./3.*X*X*X*X + 56./3.*X*X*X*Y + (44.*Y - 4.*s_pi2 + 2. - 32.*Y*Y)*X*X
                   + (-4.*Y*Y*Y - 8. - 32.*s_z3 - 80./3.*s_pi2 + 6.*Y*Y + 56./3.*s_pi2*Y)*X
                   + Y*Y*Y*Y + 6.*Y*Y*Y + (-10./3.*s_pi2 - 5.)*Y*Y + (-39. - 18.*s_pi2 + 144.*s_z3)*Y
                   + 3.*S + 187./4. - 4.*s_pi2*S + 4./45.*s_pi4 - 5.*s_pi2 - 20.*s_z3 + 48.*s_z3*S);
    const double b(-12.*X*X + (24.*Y + 24.)*X - 12.*Y*Y - 24.*Y - 12.*s_pi2);
    const double c(-80.*k.L4y + 32.*X*k.L3x + (-128.*X - 152.)*k.L3y + 152.*k.L2x*X
                   + 8.*Y*Y*k.L2y + (-16.*Y*Y - 24.)*X*X + (60.*Y*Y + (28. + 32./3.*s_pi2)*Y - 58.)*X
                   + 14./3.*Y*Y*Y*Y + 44./3.*Y*Y*Y + 8./3.*Y*Y*s_pi2 + (96.*s_z3 - 32./3.*s_pi2)*Y
                   + 32./45.*s_pi4 + 16.*s_z3 - 86./3.*s_pi2 - 2.);
    return a*k.tu + b*k.ts2 + 8.*X*X*k.tu2 + c;
  }

  double gE4(const Kin &k, double S)
  {
    const double X(k.X), Y(k.Y);
    const double a(16./9.*Y*Y*Y - 4./9.*Y*Y + 16./3.*Y*Y*S - 20.*Y + 64./9.*s_pi2*Y + 16.*Y*S
                   + 32./9.*X*X*X - 152./9.*X*X + 32./3.*X*X*S + 128./9.*s_pi2*X
                   - 2./3.*s_pi2*S + 110./9.*s_pi2 + 3401./81. - 908./27.*S - 4./9.*s_z3);
    const double c(32./9.*X*X*X - 104./9.*X*X + 32./3.*X*X*S - 152./9.*X + 32./3.*X*S
                   + 128./9.*s_pi2*X + 64./9.*s_pi2);
    return a*k.tu + c;
  }

  double gF2(const Kin &k, double S)
  { return (-92./27.*s_pi2 + 32./9.*S*S - 160./27.*S)*k.tu; }

  double gG1(const Kin &k, double)
  {
    const double X(k.X), Y(k.Y);
    const double a(14.*X*X*X*X + 28.*X*X*X + 8.*X*X*Y*Y + 56.*X*X*s_pi2 - 48.*X*X + 12.*X*X*Y
                   + 32.*X*Y*s_pi2 + 80.*s_pi2*X + 2.*Y*Y*Y*Y + 12.*Y*Y*Y - 10.*Y*Y + 8.*Y*Y*s_pi2
                   + 26.*s_pi2 + 24.*s_pi2*Y - 84.*Y + 102.);
    const double c(32.*X*X*X + 8.*X*X*Y*Y + 80.*s_pi2*X + 32.*X*Y*s_pi2 + 8.*Y*Y*X + 8.*Y*Y*Y*Y
                   + 32.*Y*Y*s_pi2 - 32.*Y*Y - 4. - 56.*Y + 24.*s_pi2);
    return a*k.tu + 8.*X*(X*X*X + X*X + 4.*s_pi2*X + 2.*s_pi2)*k.tu2
      + 2.*X*X*(X*X + 4.*s_pi2)*k.tu3 + c;
  }

  double gX4(const Kin &k, double S)
  {
    const double X(k.X), Y(k.Y);
    const double a(-32./3.*Y*s_pi2 + 32./3.*S*X*X + 16.*S*Y - 64./3.*X*s_pi2 - 16.*s_pi2
                   + 16./3.*S*Y*Y - 112./3.*S);
    const double c(32./3.*S*X*X + 32./3.*S*Y - 64./3.*X*s_pi2 - 32./3.*s_pi2);
    return a*k.tu + c;
  }

  double gX5(const Kin &k, double S)
  { return (32./9.*S*S + 32./9.*s_pi2)*k.tu; }

}

double As (double x, double S) { return Sym(gA,  x, S); }
double Bs (double x, double S) { return Sym(gB,  x, S); }
double E4s(double x, double S) { return Sym(gE4, x, S); }
double F2s(double x, double S) { return Sym(gF2, x, S); }
double G1s(double x, double S) { return Sym(gG1, x, S); }
double X4s(double x, double S) { return Sym(gX4, x, S); }
double X5s(double x, double S) { return Sym(gX5, x, S); }
double T0 (double x) { return 8.*(x/(1. - x) + (1. - x)/x); }

double Fin02(double x, double S, double sumQ4, double NF)
{
  return 2.*(sumQ4*As(x, S) + Bs(x, S) + NF*E4s(x, S) + NF*NF*F2s(x, S))/T0(x);
}

double Fin11(double x, double S, double NF)
{
  return (G1s(x, S) + NF*X4s(x, S) + NF*NF*X5s(x, S))/T0(x);
}

std::complex<double> R1(double x)
{
  // generated from derive_oneloop.py (sympy), photonic, mu^2 = s
  const double y(1. - x), X(std::log(x)), Y(std::log(y)), d(2.*(x*x + y*y));
  const double re((X*X*(1. + x*x) + X*(x*x - 4.*x + 3.) + Y*Y*(x*x - 2.*x + 2.)
                   + Y*(x*x + 2.*x) - 7.*(x*x + y*y))/d);
  const double im(s_pi*(2.*X*(1. + x*x) + 2.*Y*(x*x - 2.*x + 2.) + 2.*x*x - 2.*x + 3.)/d);
  return std::complex<double>(re, im);
}

double RF()     { return -4.; }

double F0LogFinite()
{
  // massless quark form factor, C_F^2, alpha/4pi, mu^2 = Q^2: F2 - F1^2/2 at eps^0
  const double a(-2.), b(-3.), c(-8. + s_z2), d(-16. + 1.5*s_z2 + 14./3.*s_z3),
    e(-32. + 4.*s_z2 + 7.*s_z3 + 47./20.*s_z2*s_z2);
  const double f2(1151./8. + 17./2.*s_z2 - 58.*s_z3 - 13.*s_z2*s_z2);
  return f2 - 0.5*(c*c + 2.*a*e + 2.*b*d);
}

} }
