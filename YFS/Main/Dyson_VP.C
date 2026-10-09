/*!
  \file Dyson_VP.C

  YFS: DYSON_VP - see Dyson_VP.H for the formulas and the scheme bookkeeping.
*/

#include "YFS/Main/Dyson_VP.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"
#include "MODEL/Main/Model_Base.H"

#include <cmath>
#include <sstream>

using namespace YFS;
using namespace ATOOLS;

#ifdef USING__HADALPHAQED
extern "C" {
  // MODEL/Main/hadr5x<version>.f, built with -fdefault-real-8 (see
  // MODEL/Main/Running_AlphaQED.C): E is the signed ENERGY, E < 0 spacelike
  void hadr5x_(double *e, double *st2, double *der, double *errdersta,
               double *errdersys, double *deg, double *errdegsta,
               double *errdegsys);
}
#endif

namespace {

  Settings &Main() { return Settings::GetMainSettings(); }

  //! the configurations of a massless neutral-current 2 -> 2 fermion Born
  enum class nc_config { s_same, s_flip, bhabha_same, bhabha_flip, t_only };

}


Dyson_VP::Dyson_VP(bool complete) : m_complete(complete)
{
  m_alpha0 = 1./Main()["YFS"]["1/ALPHAQED(0)"].SetDefault(137.03599976).Get<double>();
  const double amodel(MODEL::s_model->ScalarConstant("alpha_QED"));
  m_deltar = 1. - m_alpha0/amodel;
  m_mz = Flavour(kf_Z).Mass();
  m_gz = Flavour(kf_Z).Width();
  m_sw2 = MODEL::s_model->ComplexConstant("csin2_thetaW");
  m_alphas = 0.118;
#ifndef USING__HADALPHAQED
  THROW(fatal_error, "YFS: DYSON_VP needs the hadronic Delta alpha of hadr5x"
        " (SHERPA_ENABLE_ALPHAQEDHAD).");
#endif
}


Complex Dyson_VP::FermionLoop(double x, double m, double nc, double q) const
{
  /*
    (alpha0/3pi) N_c Q^2 G(x), G = 1/3 - (1 + 2m^2/x)[2 + beta ln((beta-1)/(beta+1))],
    beta = sqrt(1 - 4m^2/x), x + i0: the on-shell one-loop fermion
    vacuum polarisation (G -> ln(|x|/m^2) - 5/3 for |x| >> m^2, Im G = -pi
    (1 + 2m^2/x) beta above threshold, G(0) = 0).
  */
  if (x == 0.) return Complex(0., 0.);
  const double r(m*m/x);
  Complex bl;
  if (x < 0.) {                       // spacelike, beta > 1
    const double b(sqrt(1. - 4.*r));
    bl = Complex(b*log((b - 1.)/(b + 1.)), 0.);
  }
  else if (4.*r < 1.) {               // timelike above threshold, 0 < beta < 1
    const double b(sqrt(1. - 4.*r));
    bl = b*Complex(log((1. - b)/(1. + b)), M_PI);
  }
  else {                              // timelike below threshold, beta = i b
    const double b(sqrt(4.*r - 1.));
    bl = Complex(-2.*b*atan(1./b), 0.);
  }
  const Complex G(1./3. - (1. + 2.*r)*(2. + bl));
  return m_alpha0/(3.*M_PI)*nc*q*q*G;
}


Complex Dyson_VP::Hadronic(double x) const
{
  double re(0.);
#ifdef USING__HADALPHAQED
  double e(x < 0. ? -sqrt(-x) : sqrt(x)), st2(0.23153);
  double errsta, errsys, deg, degsta, degsys;
  hadr5x_(&e, &st2, &re, &errsta, &errsys, &deg, &degsta, &degsys);
#endif
  // Im Delta alpha_had = -(alpha0/3) R(s), perturbative R of the open quarks
  double im(0.);
  if (x > 0.) {
    const double E(sqrt(x));
    double rq(0.);
    if (E > 0.3)  rq += 3.*(4./9. + 1./9. + 1./9.);
    if (E > 3.0)  rq += 3.*4./9.;
    if (E > 10.5) rq += 3.*1./9.;
    im = -m_alpha0/3.*rq*(1. + m_alphas/M_PI);
  }
  return Complex(re, im);
}


Complex Dyson_VP::DeltaAlpha(double x) const
{
  Complex pi(Hadronic(x));
  for (kf_code kf : {kf_e, kf_mu, kf_tau})
    pi += FermionLoop(x, Flavour(kf).Mass(true), 1., 1.);
  pi += FermionLoop(x, Flavour(kf_t).Mass(true), 3., 2./3.);
  return pi;
}


Complex Dyson_VP::DysonFactor(double x) const
{
  const Complex pi(DeltaAlpha(x));
  const Complex all((1. - m_deltar)/(1. - pi) - 1.);
  return m_complete ? all : all - (pi - m_deltar);
}


Complex Dyson_VP::ZCoupling(const Flavour &f, int l) const
{
  // e^2 units: L (l = -1) (T3 - Q sw^2)/(sw cw), R (l = +1) -Q sw^2/(sw cw);
  // f the particle, so an antiparticle line passes the particle's coupling
  const Complex swcw(sqrt(m_sw2*(1. - m_sw2)));
  const double Q(f.IsAnti() ? -f.Charge() : f.Charge());
  const double T3(f.IsAnti() ? -f.IsoWeak() : f.IsoWeak());
  return (l < 0 ? T3 - Q*m_sw2 : -Q*m_sw2)/swcw;
}


Complex VP_Factors::Share(int l1, int l2, int l3, int l4) const
{
  for (const VP_Helicity &h : hel)
    if (h.l[0] == l1 && h.l[1] == l2 && h.l[2] == l3 && h.l[3] == l4) return h.r;
  return rbar;
}


VP_Factors Dyson_VP::Factors(const Flavour_Vector &fl, const Vec4D_Vector &p) const
{
  VP_Factors out;
  if (fl.size() != 4 || p.size() < 4) return out;
  const Flavour &e(fl[0]), &f(fl[2]);
  // one neutral-current pair: e- e+ -> f fbar (particle first in each pair)
  if (!e.IsLepton() || e.Charge() == 0. || e.IsAnti() || fl[1] != e.Bar() ||
      f.IsAnti() || fl[3] != f.Bar() || !f.IsFermion()) return out;
  const bool bhabha(f == e);
  const double Qe(e.Charge()), Qf(f.Charge());
  out.s = (p[0] + p[1]).Abs2();
  out.t = (p[0] - p[2]).Abs2();
  const double u((p[0] - p[3]).Abs2());
  out.Ds = DysonFactor(out.s);
  out.Dt = bhabha ? DysonFactor(out.t) : Complex(0., 0.);
  const Complex pzs(1./Complex(out.s - m_mz*m_mz, m_gz*m_mz));
  const Complex pzt(1./Complex(out.t - m_mz*m_mz, m_gz*m_mz));
  auto Fs = [&](int le, int lf, bool photon) {
    return photon ? Complex(Qe*Qf/out.s, 0.)
                  : Qe*Qf/out.s + ZCoupling(e, le)*ZCoupling(f, lf)*pzs; };
  auto Ft = [&](int le, int lep, bool photon) {
    return photon ? Complex(Qe*Qe/out.t, 0.)
                  : Qe*Qe/out.t + ZCoupling(e, le)*ZCoupling(e, lep)*pzt; };
  double norm(0.);
  for (int l1 : {-1, 1}) for (int l2 : {-1, 1}) for (int l3 : {-1, 1}) for (int l4 : {-1, 1}) {
    Complex A(0., 0.), Ag(0., 0.);
    double kin(0.);
    if (l2 == -l1 && l4 == -l3) {
      // s channel, e-(l1) e+(-l1) -> f(l3) fbar(-l3); Bhabha adds the t
      // channel when the e- keeps its helicity (both lines conserve it)
      kin = (l3 == l1 ? u : out.t);
      A  = Fs(l1, l3, false);
      Ag = out.Ds*Fs(l1, l3, true);
      if (bhabha && l3 == l1) { A += Ft(l1, l1, false); Ag += out.Dt*Ft(l1, l1, true); }
    }
    else if (bhabha && l2 == l1 && l3 == l1 && l4 == l2) {
      // t channel only: e-(l1) -> e-(l1), e+(l1) -> e+(l1)
      kin = out.s;
      A  = Ft(l1, -l1, false);
      Ag = out.Dt*Ft(l1, -l1, true);
    }
    else continue;
    VP_Helicity h;
    h.l[0] = l1; h.l[1] = l2; h.l[2] = l3; h.l[3] = l4;
    h.w = kin*kin*std::norm(A);
    h.r = std::abs(A) > 0. ? Ag/A : Complex(0., 0.);
    norm += h.w;
    out.hel.push_back(h);
  }
  if (!(norm > 0.)) return VP_Factors();
  for (VP_Helicity &h : out.hel) {
    h.w /= norm;
    out.rbar += h.w*h.r;
    out.v2   += h.w*(std::norm(1. + h.r) - 1.);
  }
  out.ok = true;
  return out;
}


std::string Dyson_VP::Describe() const
{
  std::ostringstream o;
  const Complex pz(DeltaAlpha(m_mz*m_mz)), p250(DeltaAlpha(250.*250.)),
                pt(DeltaAlpha(-10.));
  o<<(m_complete ? "complete (provider without the photon self-energy)"
                  : "resummed beyond the provider's first order")
   <<"; Delta r = 1 - alpha(0)/alpha_model = "<<m_deltar
   <<"; Pi(M_Z^2) = "<<pz<<", Pi(250^2) = "<<p250<<", Pi(-10 GeV^2) = "<<pt
   <<"; D(M_Z^2) = "<<DysonFactor(m_mz*m_mz)<<", D(250^2) = "
   <<DysonFactor(250.*250.)<<", D(-10) = "<<DysonFactor(-10.);
  return o.str();
}
