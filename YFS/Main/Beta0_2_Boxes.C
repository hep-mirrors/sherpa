/*!
  \file Beta0_2_Boxes.C

  The non-factorisable photonic two-loop beta_0^(2) of e+e- -> l+l-; see
  Beta0_2_Boxes.H.
*/

#include "YFS/Main/Beta0_2_Boxes.H"
#include "YFS/Main/Beta0_2.H"
#include "YFS/Main/Beta0_2_Boxes_Table.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_VV.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"

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
  double norm(0.);
  for (int le : {-1, 1}) for (int ll : {-1, 1}) {
    const Complex Ag(Qe*Ql/s, 0.);
    const Complex A(Ag + ZCoupling(e, le, sw2)*ZCoupling(l, ll, sw2)*pz);
    if (std::abs(A) == 0.) continue;
    const bool same(ll == le);
    const double w(sqr(same ? u : t)*std::norm(A));
    const Complex d2(a2*(same ? nsame : nopp)*Ag/A);
    f.boxhel.push_back({le, -le, ll, -ll});
    f.boxd2.push_back(d2);
    f.boxavg += w*d2.real();
    norm += w;
  }
  if (norm > 0.) f.boxavg /= norm;
  else { f.boxhel.clear(); f.boxd2.clear(); f.boxavg = 0.; }
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
  /*
    The hard function is leading power in m_e^2/|t|, m_e^2/|u|: a photon
    (nearly) along a beam of the Born point is outside it (the massified
    poles then miss YFS's by up to 1e-7 instead of 1e-11, VV_POLE_CHECK).
    YFS: BETA0_2_DIPHOTON_POWER_CUT: no term there (counted).
  */
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
