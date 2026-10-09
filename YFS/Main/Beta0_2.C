/*!
  \file Beta0_2.C

  The photonic two-loop beta_0^(2) at LL + NLL; see Beta0_2.H for the
  formulas and their derivation.
*/

#include "YFS/Main/Beta0_2.H"
#include "ATOOLS/Math/MathTools.H"

#include <cmath>
#include <sstream>

using namespace YFS;
using namespace ATOOLS;

namespace {

  const double s_zeta2(M_PI*M_PI/6.);
  const double s_zeta3(1.2020569031595942854);

  //! Z_i = theta_i Q_i with theta = +1 for the nin incoming legs
  double SignedCharge(const Flavour &fl, size_t i, size_t nin)
  { return (i < nin ? 1. : -1.)*fl.Charge(); }

  //! a crossed vertex gets C L + K only above this -t/m^2: the constant is
  //! the leading power in m^2/Q^2 (the MUonE electron line, not its muon line)
  const double s_crossed_min_q2m2(100.);

  //! one current producing or annihilating f fbar on one side
  bool SameSideVertex(const Flavour &a, const Flavour &b, double ma, double mb)
  { return a == b.Bar() && ma == mb; }

  /*
    Whether a flavour that is both incoming and outgoing is one fermion line
    through a space-like vertex (mu e -> mu e): no same-side f fbar pair and
    no flavour twice on one side, and no external photon (a line with a
    photon on it, Compton, is not one current). Bhabha (s- and t-channel
    vertices, of which only the s-channel ones are counted) and Moller (two
    possible lines per electron) are excluded.
  */
  bool UniqueCrossedLines(const Flavour_Vector &fl, size_t n, size_t nin)
  {
    for (size_t i(0); i < n; ++i) {
      if (fl[i].IsPhoton()) return false;
      if (fl[i].Charge() == 0.) continue;
      for (size_t j(i+1); j < n; ++j) {
        if (fl[j].Charge() == 0. || (i < nin) != (j < nin)) continue;
        if (fl[i] == fl[j] || fl[i] == fl[j].Bar()) return false;
      }
    }
    return true;
  }

  //! (1/2) delta_f^2 (cross) or sum_d (1/2) f_d^2 (no cross)
  Complex LeadingLog(const Photonic_Factors &f, bool cross)
  {
    if (cross) return 0.5*f.deltaf*f.deltaf;
    Complex sum(0., 0.);
    for (const Complex &fd : f.f) sum += 0.5*fd*fd;
    return sum;
  }

  //! exp(delta_f)(1 + NLL) - 1 - delta_f, per pair without the cross products
  Complex Exponentiated(const Photonic_Factors &f, bool cross)
  {
    if (cross) return exp(f.deltaf)*(1. + f.nll) - 1. - f.deltaf;
    Complex sum(0., 0.);
    for (size_t d(0); d < f.f.size(); ++d)
      sum += exp(f.f[d])*(1. + f.n[d]) - 1. - f.f[d];
    return sum;
  }

  //! the factor w multiplies in the delta_f w product
  Complex WeakProductFactor(const Photonic_Factors &f, const Beta02_Config &cfg)
  {
    if (cfg.mode == beta02::exp && cfg.cross)
      return exp(f.deltaf)*(1. + f.nll) - 1.;
    return f.deltaf;
  }

  //! delta^(2) for a delta_f and a w of the same kind (both complex or both real)
  Complex Assemble(const Photonic_Factors &f, const Complex &w,
                   const Beta02_Config &cfg)
  {
    Complex d2(0., 0.);
    switch (cfg.mode) {
    case beta02::none:  return d2;
    case beta02::ll:    d2 = LeadingLog(f, cfg.cross); break;
    case beta02::llnll: d2 = LeadingLog(f, cfg.cross) + f.nll; break;
    case beta02::two_loop: d2 = LeadingLog(f, cfg.cross) + f.nll + f.cst; break;
    case beta02::exp:   d2 = Exponentiated(f, cfg.cross); break;
    }
    if (cfg.weak) d2 += w*WeakProductFactor(f, cfg);
    return d2;
  }

  //! the same factors with every imaginary part dropped
  Photonic_Factors RealPart(const Photonic_Factors &f)
  {
    Photonic_Factors r(f);
    r.deltaf = Complex(f.deltaf.real(), 0.);
    r.nll    = Complex(f.nll.real(), 0.);
    for (Complex &x : r.f) x = Complex(x.real(), 0.);
    for (Complex &x : r.n) x = Complex(x.real(), 0.);
    r.cst = Complex(f.cst.real(), 0.);
    return r;
  }

}


double YFS::Beta02NLLCoefficient()
{
  return 3./32. - 0.75*s_zeta2 + 1.5*s_zeta3;
}


double YFS::Beta02TwoLoopConstant()
{
  // Penin hep-ph/0508127 eq. (21), the constant of the massive form factor
  const double pi2(M_PI*M_PI);
  const double K(11./8. + 17./32.*pi2 - 2.25*s_zeta3 - 2./45.*pi2*pi2
                 - 0.5*pi2*log(2.));
  return Beta02NLLCoefficient() + K;
}


Complex YFS::Beta02CnuA(double svar, double m1, double m2)
{
  const double mm(m1*m2);
  const Complex Mas12(mm, 0.), Nu((-svar + m1*m1 + m2*m2)/2., 0.);
  const Complex xlam(sqrt((Nu - Mas12)*(Nu + Mas12)));
  const Complex z(svar > 0. ? Mas12/(Nu - xlam) : (Nu + xlam)/Mas12);
  // KKbvir::CDLN(z, +1): +i pi on the negative real axis
  const Complex lz((z.imag() == 0. && z.real() <= 0.)
                   ? log(-z) + Complex(0., M_PI) : log(z));
  return Nu/xlam * lz;
}


Photonic_Factors YFS::PhotonicFactors(const Flavour_Vector &fl,
                                      const Vec4D_Vector &p,
                                      const Beta02_Config &cfg, size_t nin)
{
  Photonic_Factors out;
  const size_t n(Min(p.size(), fl.size()));
  const double a2pi(cfg.alpha/(2.*M_PI)), api(cfg.alpha/M_PI);
  const double cnll(Beta02NLLCoefficient());
  // computed for every mode: the factors are shared by the two columns
  const double c2(Beta02TwoLoopConstant());
  const bool crossedlines(UniqueCrossedLines(fl, n, nin));
  bool any(false);
  for (size_t i(0); i < n; ++i) {
    const double Qi(fl[i].Charge());
    if (Qi == 0.) continue;
    const double mi(fl[i].Mass());
    if (!(mi > 0.)) return Photonic_Factors();
    for (size_t j(i+1); j < n; ++j) {
      const double Qj(fl[j].Charge());
      if (Qj == 0.) continue;
      const double mj(fl[j].Mass());
      if (!(mj > 0.)) return Photonic_Factors();
      const bool same((i < nin) == (j < nin));
      if (cfg.pairs == beta02pairs::same_side && !same) continue;
      const double eta(-SignedCharge(fl[i], i, nin)*SignedCharge(fl[j], j, nin));
      const double sij(same ? (p[i]+p[j]).Abs2() : (p[i]-p[j]).Abs2());
      const Complex K(a2pi*eta*(Beta02CnuA(sij, mi, mj) - 1.));
      // the phase convention of B_V: Im delta shifts by -Im K (c - 1/2)
      const Complex fij(K - Complex(0., K.imag()*(cfg.imc - 0.5)));
      const Complex nll(api*cnll*((Qi*Qi + Qj*Qj)*K
                                  + a2pi*eta*(Qi*Qi - Qj*Qj)*log(mj/mi)));
      // a vertex: one current, a same-side f fbar pair or, across, one line
      // in and out with Q^2 = -s_ij (Penin's Euclidean case, no continuation)
      const bool vertex(same ? SameSideVertex(fl[i], fl[j], mi, mj)
                        : crossedlines && fl[i] == fl[j]
                          && -sij > s_crossed_min_q2m2*mi*mi);
      const Complex cij(vertex ? api*api*Qi*Qi*Qj*Qj*c2 : 0.);
      out.f.push_back(fij);
      out.n.push_back(nll);
      out.c.push_back(cij);
      out.deltaf += fij;
      out.nll    += nll;
      out.cst    += cij;
      any = true;
    }
  }
  out.ok = any;
  return out;
}


Complex YFS::Beta02Delta(const Photonic_Factors &f, const Complex &delta,
                         const Beta02_Config &cfg)
{
  if (!f.ok || !cfg.On()) return Complex(0., 0.);
  return Assemble(f, delta - f.deltaf, cfg);
}


double YFS::Beta02DeltaReal(const Photonic_Factors &f, double delta,
                            const Beta02_Config &cfg)
{
  if (!f.ok || !cfg.On()) return 0.;
  const Photonic_Factors r(RealPart(f));
  // without helicities the two-loop boxes enter Born-weighted
  return Assemble(r, Complex(delta - r.deltaf.real(), 0.), cfg).real()
    + (cfg.mode == beta02::two_loop ? f.boxavg : 0.);
}


Complex YFS::Beta02VertexFactor(const Complex &F1, double alpha, double Q2,
                                beta02::code mode)
{
  const Complex nll(2.*Beta02NLLCoefficient()*(alpha/M_PI)*Q2*F1);
  switch (mode) {
  case beta02::none:  return F1;
  case beta02::ll:    return F1 + 0.5*F1*F1;
  case beta02::llnll: return F1 + 0.5*F1*F1 + nll;
  case beta02::exp:   return exp(F1)*(1. + nll) - 1.;
  case beta02::two_loop:
    return F1 + 0.5*F1*F1 + nll
      + sqr(alpha/M_PI)*Q2*Q2*Beta02TwoLoopConstant();
  }
  return F1;
}


std::string YFS::Beta02Describe(const Beta02_Config &cfg)
{
  std::ostringstream o;
  o<<"mode "<<cfg.mode<<", pairs "<<cfg.pairs
   <<", cross products "<<(cfg.cross ? "on" : "off")
   <<", delta_f w "<<(cfg.weak ? "on" : "off")
   <<", 1/alpha "<<1./cfg.alpha<<", c "<<cfg.imc
   <<", c_NLL "<<Beta02NLLCoefficient();
  return o.str();
}
