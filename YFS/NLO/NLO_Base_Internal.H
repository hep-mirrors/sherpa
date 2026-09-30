#ifndef YFS_NLO_NLO_Base_Internal_H
#define YFS_NLO_NLO_Base_Internal_H

// Private to the NLO_*.C files of YFS/NLO (not installed): the helpers
// that more than one of them uses, moved out of NLO_Base.C unchanged.
// They are inline, not static, so that each function-local static (the
// cached YFS settings) exists once in the library rather than once per
// translation unit. The using-directives are those of the NLO_*.C files,
// which the helper bodies were written against.

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "YFS/Main/Photon.H"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

/*
  YFS: RV_MODE - how the squared-level real-virtual (NLO_Part with E) is built.

  0 (default): the legacy RV, CalculateRealVirtual(Vec4D): per photon
    [r_RV flux - (B~ rtree flux + S~ B v)/kappa]/S~_crude, the whole
    beta_1^(2) - beta_1^(1) = R v_{n+1} - S~ B v, with its own point, crude,
    flux and coupling, and with E present REAL_COMBINE's product and
    VIRTUAL_COMBINE switched off. IR unsafe: the (n+1) loop runs at
    mu = 1 GeV against the Born's mu^2 = s, and the full-EW finite part is
    mu dependent beyond its IR pole (+0.4% G_mu, +0.8% alpha(0) per event
    whose photons are all soft).

  1: the remainder. With v = V_sub/B of the event, v_{n+1} the same object at
    photon j's (n+1)-body point and rho_j = R_j/(crude_j B) the photon's
    exact-over-crude ratio, the exact one-photon O(alpha^2) weight
      beta_0 (1 + v) + [R (1 + v_{n+1}) - S B (1 + v)]/crude
    is identically
      (1 + v)(1 + delta) + rho (v_{n+1} - v),
    delta = (R - S B)/(crude B) the real's bracket.
*/
/*
  YFS: RR_MODE - how the squared-level double real (NLO_Part with W) is built.
  0 (default): the legacy assembly (CalculateRealReal(k1,k2), steered by
    RR_CONVENTIONS), which also switches REAL_COMBINE's product off.
  1: exact at O(alpha^2) for every photon pair. f_k is photon k's factor in
    REAL_COMBINE's product (its bracket plus its IF ratio, = kappa r F/(c B),
    the exact over crude at its own single-real point), and
        rho_ij = kappa^2 R_2(P) F_ij / (c_i c_j B)
    the pair's exact over crude at its point P (EventPairPoint: the event,
    the other photons reduced), c the singles' densities carried to P by the
    ratio of crude eikonals, F_ij the pair's flux (P_beams - K_ISR)^2/P_beams^2.
    The pair adds
        RR_ij = B (rho_ij - f_i f_j) prod_{k != i,j} f_k,
    so that with two photons product + RR = rho_ij, the exact two-photon ME
    over the generator's density, event by event. RR does NOT vanish when
    one photon goes soft next to a hard one: the exact R_2 attaches the soft
    photon to the legs after the hard emission, the product (and the form
    factor) to the generation legs. That soft x hard log is cancelled by the
    real-virtual of the hard photon (RV_MODE 1, CalculateRealVirtualRemainder,
    the form-factor difference at the RR's cut), so RR_MODE 1 needs RV_MODE 1
    and RR_SOFT_CUT > 0, and the test is sigma independent of RR_SOFT_CUT and
    IR_CUTOFF, not an event-level soft limit (notes sec. 11.12).
    Settings: RR_Generator OpenLoops (Comix's R_2 is unreliable for soft
    photons), RR_SOFT_CUT 1e-4 (OpenLoops' R_2 for two ultra-soft collinear
    initial-state photons, x ~ 1e-7, depends on the frame). Born photons'
    channel weights (REAL_BORN_PHOTON_MULTICHANNEL) are not applied to R_2
    (e+e- -> gamma gamma needs them). RR_PROBE 1 prints "@@@ RRPAIR" per pair
    (2 adds P and the Born point) and "@@@ RRID" per two-photon event.
*/
inline int RRMode() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_MODE"]
                     .SetDefault(0).Get<int>());
  return m;
}

inline int RVMode() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_MODE"]
                     .SetDefault(0).Get<int>());
  return m;
}

// YFS: RV_LOOP_FRAME, documented at RVLoopPoint
inline int RVLoopFrame() {
  static const int m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_LOOP_FRAME"]
                     .SetDefault(1).Get<int>());
  return m;
}

/*
  The frame the one-loop provider is called in (RV_MODE 1 always, the Born's
  virtual with YFS: VIRTUAL_CANONICAL_FRAME). OpenLoops' one-loop amplitudes
  are not Lorentz invariant to better than ~1% when a leg is almost, but not
  exactly, on a coordinate axis - measured with pyol at Sherpa's own points

  So: boost to the rest frame of the two beams, put the beams exactly back to
  back on the z axis, then turn the whole point by
  a fixed generic rotation, so that no leg lies near a coordinate axis. The
  squared amplitudes are invariant, so nothing else changes.
*/
inline Vec4D_Vector BeamRestFrameOnZ(const Vec4D_Vector &p)
{
  Vec4D_Vector q(p);
  if (q.size() < 2) return q;
  Poincare cms(q[0] + q[1]);
  for (Vec4D &v : q) cms.Boost(v);
  if (sqr(q[0][1]) + sqr(q[0][2]) > 0. || q[0][3] < 0.) {
    const Poincare rot(q[0], Vec4D(0., 0., 0., 1.));
    Vec4D a(q[0]), b(q[0]);
    rot.Rotate(a);
    rot.RotateBack(b);
    const bool fwd(sqr(a[1]) + sqr(a[2]) <= sqr(b[1]) + sqr(b[2]));
    for (Vec4D &v : q) { if (fwd) rot.Rotate(v); else rot.RotateBack(v); }
  }
  q[0] = Vec4D(q[0][0], 0., 0.,  Vec3D(q[0]).Abs());
  q[1] = Vec4D(q[1][0], 0., 0., -Vec3D(q[1]).Abs());
  return q;
}

inline Vec4D_Vector CanonicalBeamFrame(const Vec4D_Vector &p)
{
  Vec4D_Vector q(BeamRestFrameOnZ(p));
  // fixed generic rotation (z-y-z Euler angles 0.3, 0.7, 1.1)
  static const double ca(cos(0.3)), sa(sin(0.3)), cb(cos(0.7)), sb(sin(0.7)),
                      cc(cos(1.1)), sc(sin(1.1));
  for (Vec4D &v : q) {
    double x(v[1]), y(v[2]), z(v[3]);
    double x1(ca*x - sa*y), y1(sa*x + ca*y);
    double x2(cb*x1 + sb*z), z2(-sb*x1 + cb*z);
    v = Vec4D(v[0], cc*x2 - sc*y1, sc*x2 + cc*y1, z2);
  }
  return q;
}

template <typename Get>
inline void HardestBetas(const YFS::Photon_Vector &ph, Get get,
                         double &h1, double &h2) {
  std::vector<const YFS::Photon *> s;
  s.reserve(ph.size());
  for (const YFS::Photon &k : ph) s.push_back(&k);
  std::sort(s.begin(), s.end(),
            [](const YFS::Photon *a, const YFS::Photon *b) {
              return a->E() > b->E();
            });
  h1 = s.empty() ? 0. : get(*s[0]);
  h2 = h1 + (s.size() < 2 ? 0. : get(*s[1]));
}

#endif
