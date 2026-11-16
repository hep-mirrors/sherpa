/*!
  \file Beta0_2_Boxes.C

  The non-factorisable photonic two-loop beta_0^(2) of e+e- -> l+l-; see
  Beta0_2_Boxes.H.
*/

#include "YFS/Main/Beta0_2_Boxes.H"
#include "YFS/Main/Beta0_2.H"
#include "YFS/Main/Beta0_2_Boxes_Table.H"
#include "YFS/Main/Beta0_2_Boxes_Bhabha_Table.H"
#include "YFS/Main/Beta0_2_Boxes_MuE_Table.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_VV.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <cmath>

using namespace YFS;
using namespace ATOOLS;

namespace {

  //! Catmull-Rom interpolation of one table column at y (clamped to the grid)
  double Interpolate(const double *v, double y)
  {
    using namespace B02Table;
    const double u(Max(0., Min(double(s_n - 1), (y - s_y0)/s_dy)));
    const int i(Min(s_n - 2, int(u)));
    const double t(u - i);
    const double p0(v[Max(0, i - 1)]), p1(v[i]), p2(v[i + 1]),
                 p3(v[Min(s_n - 1, i + 2)]);
    return p1 + 0.5*t*(p2 - p0 + t*(2.*p0 - 5.*p1 + 4.*p2 - p3
                                    + t*(3.*(p1 - p2) + p3 - p0)));
  }

  //! the Z coupling in e^2 units, L (l = -1) or R (l = +1), as Dyson_VP
  Complex ZCoupling(const Flavour &f, int l, const Complex &sw2)
  {
    const Complex swcw(sqrt(sw2*(1. - sw2)));
    return (l < 0 ? f.IsoWeak() - f.Charge()*sw2 : -f.Charge()*sw2)/swcw;
  }

  //! Bhabha boxes only at leading power: min(|t|, |u|) above this times m_e^2
  const double s_bhabha_min_q2m2(100.);
  //! l e -> l e remainder only at leading power: |t| above this times m_e^2
  const double s_mue_min_q2m2(100.);

  //! the three Bhabha helicity classes of Beta0_2_Boxes_Bhabha_Table.H
  enum class BhabhaClass { s_only, s_and_t, t_only };

  //! N(class, x) of Bhabha in (alpha/2pi)^2
  Complex TwoLoopPhotonicNBhabha(BhabhaClass c, double x)
  {
    using namespace B02BhabhaTable;
    if (!(x > 0. && x < 1.)) return Complex(0., 0.);
    if (c == BhabhaClass::s_only) return TwoLoopPhotonicN(false, x);
    const double y(log(x/(1. - x)));
    return c == BhabhaClass::s_and_t
      ? Complex(Interpolate(s_st_re, y), Interpolate(s_st_im, y))
      : Complex(Interpolate(s_t_re, y), Interpolate(s_t_im, y));
  }

  //! one helicity configuration's photon and full reduced amplitudes, Born weight
  struct Box_Term {
    std::vector<int> hel;
    BhabhaClass cls;
    Complex Ag, A;
    double spinor2;
  };

  std::vector<Box_Term> BhabhaTerms(const Flavour &e, double s, double t,
                                    double u, const Complex &ps,
                                    const Complex &pt, const Complex &sw2)
  {
    const double Q2(sqr(e.Charge()));
    std::vector<Box_Term> terms;
    for (int le : {-1, 1}) {
      const Complex gg(ZCoupling(e, le, sw2)*ZCoupling(e, le, sw2));
      const Complex go(ZCoupling(e, le, sw2)*ZCoupling(e, -le, sw2));
      const Complex As(Q2/s, 0.), Ast(Q2*(1./s + 1./t), 0.), At(Q2/t, 0.);
      terms.push_back({{le, -le, -le, le}, BhabhaClass::s_only, As, As + go*ps, t*t});
      terms.push_back({{le, -le, le, -le}, BhabhaClass::s_and_t, Ast,
                       Ast + gg*(ps + pt), u*u});
      terms.push_back({{le, le, le, le}, BhabhaClass::t_only, At, At + go*pt, s*s});
    }
    return terms;
  }

  /*
    The NLL terms of the pairs with one leg in and one out. N holds them for the
    photon exchange (it is everything beyond the two s-channel form factors), so
    each helicity's box term is (a2 N - this) A_gamma/A, and they stay on the Z share.
  */
  Complex OppositeSideNLL(const Photonic_Factors &f)
  {
    Complex sum(0., 0.);
    for (size_t d(0); d < f.n.size() && d < f.opposite.size(); ++d)
      if (f.opposite[d]) sum += f.n[d];
    return sum;
  }

  //! the Bhabha boxes of AddPhotonicTwoLoopBoxes
  void AddBhabhaBoxes(Photonic_Factors &f, const Flavour &e, const Vec4D_Vector &p,
                      double alpha)
  {
    const double s((p[0] + p[1]).Abs2()), t((p[0] - p[2]).Abs2()),
                 u((p[0] - p[3]).Abs2());
    if (!(s > 0.)) return;
    if (Min(-t, -u) < s_bhabha_min_q2m2*sqr(e.Mass())) return;
    const double x(Max(0., Min(1., -t/s)));
    const double mz(Flavour(kf_Z).Mass()), gz(Flavour(kf_Z).Width());
    const Complex sw2(MODEL::s_model->ComplexConstant("csin2_thetaW"));
    const Complex ps(1./Complex(s - mz*mz, gz*mz)), pt(1./Complex(t - mz*mz, gz*mz));
    const double a2(sqr(alpha/(2.*M_PI)));
    const Complex nllx(OppositeSideNLL(f));
    double norm(0.);
    for (const Box_Term &b : BhabhaTerms(e, s, t, u, ps, pt, sw2)) {
      if (std::abs(b.A) == 0.) continue;
      const double w(b.spinor2*std::norm(b.A));
      const Complex d2((a2*TwoLoopPhotonicNBhabha(b.cls, x) - nllx)*b.Ag/b.A);
      f.boxhel.push_back(b.hel);
      f.boxd2.push_back(d2);
      f.boxavg += w*d2.real();
      norm += w;
    }
    if (norm > 0.) f.boxavg /= norm;
    else { f.boxhel.clear(); f.boxd2.clear(); f.boxavg = 0.; }
  }

}


Complex YFS::TwoLoopPhotonicN(bool same, double x)
{
  using namespace B02Table;
  if (!(x > 0. && x < 1.)) return Complex(0., 0.);
  const double y(log(x/(1. - x)));
  return same ? Complex(Interpolate(s_same_re, y), Interpolate(s_same_im, y))
              : Complex(Interpolate(s_opp_re, y), Interpolate(s_opp_im, y));
}


void YFS::AddPhotonicTwoLoopBoxes(Photonic_Factors &f, const Flavour_Vector &fl,
                                  const Vec4D_Vector &p, double alpha)
{
  f.boxhel.clear(); f.boxd2.clear(); f.boxavg = 0.;
  if (fl.size() != 4 || p.size() < 4) return;
  const Flavour &e(fl[0]), &l(fl[2]);
  if (e.Kfcode() != kf_e || e.IsAnti() || fl[1] != e.Bar()) return;
  if (l == e && fl[3] == e.Bar()) { AddBhabhaBoxes(f, e, p, alpha); return; }
  if ((l.Kfcode() != kf_mu && l.Kfcode() != kf_tau) || l.IsAnti() ||
      fl[3] != l.Bar()) return;
  const double s((p[0] + p[1]).Abs2()), t((p[0] - p[2]).Abs2()),
               u((p[0] - p[3]).Abs2());
  if (!(s > 0.)) return;
  const double x(Max(0., Min(1., -t/s)));
  const Complex nsame(TwoLoopPhotonicN(true, x)), nopp(TwoLoopPhotonicN(false, x));
  const double mz(Flavour(kf_Z).Mass()), gz(Flavour(kf_Z).Width());
  const Complex sw2(MODEL::s_model->ComplexConstant("csin2_thetaW"));
  const Complex pz(1./Complex(s - mz*mz, gz*mz));
  const double Qe(e.Charge()), Ql(l.Charge());
  const double a2(sqr(alpha/(2.*M_PI)));
  const Complex nllx(OppositeSideNLL(f));
  double norm(0.);
  for (int le : {-1, 1}) for (int ll : {-1, 1}) {
    const Complex Ag(Qe*Ql/s, 0.);
    const Complex A(Ag + ZCoupling(e, le, sw2)*ZCoupling(l, ll, sw2)*pz);
    if (std::abs(A) == 0.) continue;
    const bool same(ll == le);
    const double w(sqr(same ? u : t)*std::norm(A));
    const Complex d2((a2*(same ? nsame : nopp) - nllx)*Ag/A);
    f.boxhel.push_back({le, -le, ll, -ll});
    f.boxd2.push_back(d2);
    f.boxavg += w*d2.real();
    norm += w;
  }
  if (norm > 0.) f.boxavg /= norm;
  else { f.boxhel.clear(); f.boxd2.clear(); f.boxavg = 0.; }
}


namespace {

  //! Catmull-Rom interpolation of v[0 .. n-1] on the uniform grid y0 + i dy
  double CatmullRom(const double *v, int n, double y0, double dy, double y)
  {
    const double u(Max(0., Min(double(n - 1), (y - y0)/dy)));
    const int i(Min(n - 2, int(u)));
    const double t(u - i);
    const double p0(v[Max(0, i - 1)]), p1(v[i]), p2(v[i + 1]), p3(v[Min(n - 1, i + 2)]);
    return p1 + 0.5*t*(p2 - p0 + t*(2.*p0 - 5.*p1 + 4.*p2 - p3
                                    + t*(3.*(p1 - p2) + p3 - p0)));
  }

  //! one column of Beta0_2_Boxes_MuE_Table.H at (s, nu), s = s/M^2: Catmull-Rom
  //! in nu, linear between the s nodes, both clamped to the grid
  double MuETable(const double (*col)[B02MuETable::s_nnu], double s, double nu)
  {
    using namespace B02MuETable;
    const double sc(Max(s_s[0], Min(s_s[s_ns - 1], s)));
    int j(0);
    while (j < s_ns - 2 && sc > s_s[j + 1]) ++j;
    const double w((sc - s_s[j])/(s_s[j + 1] - s_s[j]));
    return (1. - w)*CatmullRom(col[j], s_nnu, s_nu0, s_dnu, nu)
      + w*CatmullRom(col[j + 1], s_nnu, s_nu0, s_dnu, nu);
  }

  //! the leg indices of l e -> l e: {e in, l in, e out, l out}, or empty
  std::vector<size_t> MuELegs(const Flavour_Vector &fl)
  {
    if (fl.size() != 4) return {};
    auto light = [](const Flavour &f) { return f.Kfcode() == kf_e; };
    auto heavy = [](const Flavour &f) { return f.Kfcode() == kf_mu || f.Kfcode() == kf_tau; };
    for (size_t ei : {0, 1}) for (size_t eo : {2, 3}) {
      const size_t li(1 - ei), lo(5 - eo);
      if (light(fl[ei]) && fl[eo] == fl[ei] && heavy(fl[li]) && fl[lo] == fl[li])
        return {ei, li, eo, lo};
    }
    return {};
  }

}


double YFS::MuEPhotonicTwoLoopD(double s, double t, double qQ)
{
  if (!(s > 1.) || !(t < 0.)) return 0.;
  const double tmin(-sqr(s - 1.)/s);
  const double xi(Max(1.e-12, Min(1. - 1.e-12, t/tmin)));
  const double nu(log(xi/(1. - xi)));
  return MuETable(B02MuETable::s_even, s, nu) + qQ*MuETable(B02MuETable::s_odd, s, nu);
}


YFS::MuEBox YFS::AddMuEPhotonicTwoLoop(Photonic_Factors &f, const Flavour_Vector &fl,
                                      const Vec4D_Vector &p, double alpha)
{
  using namespace B02MuETable;
  f.muebox = 0.;
  f.mueboxapplied = false;
  const std::vector<size_t> leg(MuELegs(fl));
  if (leg.empty() || p.size() != 4) return MuEBox::not_mue;
  const double M2(sqr(fl[leg[1]].Mass()));
  if (!(M2 > 0.)) return MuEBox::not_mue;
  // s, t in units of the heavy mass, at m_e = 0
  const double s((p[leg[0]] + p[leg[1]]).Abs2()/M2), t((p[leg[0]] - p[leg[2]]).Abs2()/M2);
  if (!(s >= s_s[0] && s <= s_s[s_ns - 1])) return MuEBox::outside_s;
  if (!(-t*M2 > s_mue_min_q2m2*sqr(fl[leg[0]].Mass()))) return MuEBox::below_power;
  const double xi(t/(-sqr(s - 1.)/s));
  if (!(xi > 0.) || log(xi/Max(1.e-300, 1. - xi)) < s_nu0) return MuEBox::below_power;
  // physical charges: Q_e Q_l = -1 for mu+ e-, +1 for mu- e-
  const double qQ(fl[leg[0]].Charge()*fl[leg[1]].Charge());
  f.muebox = sqr(alpha/M_PI)*MuEPhotonicTwoLoopD(s, t, qQ);
  f.boxavg += f.muebox;
  f.mueboxapplied = true;
  static bool cited(false);
  if (!cited) {
    cited = true;
    rpa->gen.AddCitation(1, "The photonic two-loop l e -> l e amplitude (YFS: BETA0_2_MUE_BOXES) is"
                         " tabulated from \\cite{Broggio:2022htr}, based on"
                         " \\cite{Bonciani:2021okt,Mastrolia:2017pfy,DiVita:2018nnh}.");
  }
  return MuEBox::applied;
}


double YFS::MuEBoxHelicityShare(const Photonic_Factors &f, const std::vector<Complex> &a0,
                                const std::vector<Complex> &delta)
{
  if (f.muebox == 0. || a0.empty() || a0.size() != delta.size()) return 0.;
  double s0(0.), s1(0.), s2(0.);
  for (size_t c(0); c < a0.size(); ++c) {
    const double w(std::norm(a0[c])), r(delta[c].real());
    s0 += w; s1 += w*r; s2 += w*r*r;
  }
  if (!(s0 > 0.)) return 0.;
  return f.muebox - (s2/s0 - sqr(s1/s0));
}


Complex YFS::Beta02BoxFor(const Photonic_Factors &f, const std::vector<int> &l)
{
  for (size_t c(0); c < f.boxhel.size(); ++c)
    if (f.boxhel[c] == l) return f.boxd2[c];
  return Complex(0., 0.);
}


/*
  e+e- -> gamma gamma: the AGT hard function (EXTRA_XS/NNLO) in BETA0_2's
  scheme. PhotonicFactors treats the e+e- pair as an s-channel vertex, which
  is right for the LL, the per-leg NLL (C_NLL = -H/8) and the massification
  constant, and wrong only for the hard function; the difference is
    Delta(x) = Fin02(x) - f0/2 - 2 r_F Re r1(x) + r_F^2 ,   mu^2 = s,
  0 for a plain vertex (r1 -> r_F, Fin02 -> f0/2 + r_F^2). The one loop
  squared (Fin11) is ACRAIC's own |delta|^2 per helicity and is not added.
  Derivation and checks: yfs-helpers/diphoton-two-loop.
*/
double YFS::AddPhotonicTwoLoopDiphoton(Photonic_Factors &f,
                                       PHASIC::Virtual_ME2_Base *vv,
                                       const Vec4D_Vector &p, double alpha,
                                       double me, double powercut,
                                       bool &provided, bool &skipped)
{
  provided = skipped = false;
  EXTRAXS::ee_aa_QED_VV *aa(dynamic_cast<EXTRAXS::ee_aa_QED_VV*>(vv));
  if (aa == NULL) return 0.;
  provided = true;
  if (p.size() < 4) return 0.;

  const double t(std::abs((p[0] - p[2]).Abs2())), u(std::abs((p[0] - p[3]).Abs2()));
  if (powercut > 0. && me > 0. && me*me > powercut*Min(t, u)) { skipped = true; return 0.; }
  aa->Calc(p);
  if (!aa->Valid()) return 0.;
  namespace AGT = EXTRAXS::ee_aa_QED_TwoLoop;
  const double rf(AGT::RF());
  const double delta(aa->Fin02() - 0.5*AGT::F0LogFinite()
                     - 2.*rf*aa->R1().real() + rf*rf);
  const double d(0.5*sqr(alpha/(2.*M_PI))*delta);
  f.cst += Complex(d, 0.);
  return d;
}
