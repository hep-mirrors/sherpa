#ifndef YFS_NLO_NLO_Base_Internal_H
#define YFS_NLO_NLO_Base_Internal_H

// Private to the NLO_*.C files (not installed): helpers used by more than one
// of them. Inline so that each function-local static (a cached YFS setting)
// exists once per library.

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
  With v = V_sub/B of the event, v_{n+1} the same at photon j's (n+1)-body
  point and rho_j = R_j/(crude_j B), the one-photon O(alpha^2) weight
    beta_0 (1 + v) + [R (1 + v_{n+1}) - S B (1 + v)]/crude
  equals (1 + v)(1 + delta) + rho (v_{n+1} - v), delta = (R - S B)/(crude B).
  legacy (0): CalculateRealVirtual(Vec4D), with its own point, crude, flux and
    coupling; switches REAL_COMBINE's product off. IR unsafe (the loop's mu
    differs from the Born's).
  remainder (1): v_{n+1} - v = lt - V/B, YFS subtraction on the generation legs.
  remainder_own_legs (2): v_{n+1} = lt - B(point legs), plus the exponent's leg
    difference Y_FF+IF(point legs) - Y_FF+IF(generation legs). Differs from
    mode 1 only for final-state photons. ACRAIC's averaged real-virtual uses
    this convention (NLO_Base::m_rvdv).
*/
/*
  YFS: RR_MODE - how the squared-level double real (NLO_Part with W) is built.
  f_k is photon k's factor in REAL_COMBINE's product, rho_ij = kappa^2 R_2(P)
  F_ij/(c_i c_j B) the pair's exact over crude at its point P (EventPairPoint),
  F_ij the pair's flux.
  legacy (0): CalculateRealReal(k1,k2), steered by RR_CONVENTIONS; switches
    REAL_COMBINE's product off.
  exact (1): adds RR_ij = B (rho_ij - f_i f_j) prod_{k != i,j} f_k. Needs RV_MODE
    1, RR_SOFT_CUT > 0 and RR_Generator OpenLoops; not for IFI_Real false.
  singles_product (2): adds RR_ij = B (rho_ij - rho_i rho_j) prod_{k != i,j} f_k,
    rho_k = RRPhotonInfo::rho. Right in both IFI_Real schemes; taken over every
    generated pair (RR_SOFT_CUT 0).
  v x RR and v x R_n are O(alpha^3) and left out of the squared-level weight
  (YFS_Handler::CalculateNLO).
*/
inline rrmode::code RRMode() {
  static const rrmode::code m(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_MODE"]
                              .SetDefault(rrmode::singles_product).Get<rrmode::code>());
  return m;
}

//! RR_MODE exact and singles_product: the remainder is built from the
//! single reals' RRPhotonInfo and sits on top of REAL_COMBINE's product
inline bool RRFromSingles() {
  return RRMode() == rrmode::exact || RRMode() == rrmode::singles_product;
}

// YFS: REAL_COMBINE, documented where CalculateNLO's real uses it
inline realcombine::code RealCombineMode() {
  static const realcombine::code m(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_COMBINE"]
                                   .SetDefault(realcombine::product).Get<realcombine::code>());
  return m;
}

inline rvmode::code RVMode() {
  static const rvmode::code m(ATOOLS::Settings::GetMainSettings()["YFS"]["RV_MODE"]
                              .SetDefault(rvmode::remainder).Get<rvmode::code>());
  return m;
}

//! RV_MODE remainder and remainder_own_legs: the real-virtual as the
//! remainder beyond (1 + v) x real, on the real's own point
inline bool RVIsRemainder() {
  return RVMode() == rvmode::remainder || RVMode() == rvmode::remainder_own_legs;
}

// YFS: RV_LOOP_FRAME, documented at RVLoopPoint
inline rvloopframe::code RVLoopFrame() {
  static const rvloopframe::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
    ["RV_LOOP_FRAME"].SetDefault(rvloopframe::tilted).Get<rvloopframe::code>());
  return m;
}

// YFS: REAL_SUB_EIK, documented where it is used in NLO_Base::RealSubtraction
// and the SubEik* helpers (NLO_Real.C)
inline realsubeik::code RealSubEik() {
  static const realsubeik::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
    ["REAL_SUB_EIK"].SetDefault(realsubeik::multichannel_prefsr).Get<realsubeik::code>());
  return m;
}

/*
  YFS: REAL_FSR_FLUX, the flux a final-state photon's real carries (name or
  old integer; NLO_Base::RealFlux):
    event_flux (0): CalculateFlux's, as for an initial-state photon;
    no_flux (1, default): none;
    flux_squared (2): the event flux squared, (m_ff^2/s')^2;
    on_crude (3): no flux on r; F = (q + k)^2/q^2 on the crude's final-state
      part instead, q the post-emission pair;
    own_dipole (4): the flux of the photon's own pair, (Q_D - k)^2/Q_D^2.
*/
inline realfsrflux::code RealFSRFlux() {
  static const realfsrflux::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
    ["REAL_FSR_FLUX"].SetDefault(realfsrflux::no_flux).Get<realfsrflux::code>());
  return m;
}

// YFS: REAL_FSR_MAP, documented in YFS_Base::RegisterDefaults
inline realfsrmap::code RealFSRMap() {
  static const realfsrmap::code m(ATOOLS::Settings::GetMainSettings()["YFS"]
                                  ["REAL_FSR_MAP"].Get<realfsrmap::code>());
  return m;
}

// YFS: REAL_ALPHA0 (default true), documented in NLO_Base::CalculateReal:
// beta_1 at alpha(0) rather than at the model's alpha
inline bool RealAlpha0() {
  static const bool m(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_ALPHA0"]
                      .SetDefault(true).Get<bool>());
  return m;
}

/*
  The frame the one-loop provider is called in (RV_MODE 1 always, the Born's
  virtual with YFS: VIRTUAL_CANONICAL_FRAME). Loop amplitudes are not Lorentz
  invariant at the percent level for a leg near a coordinate axis, so: boost to
  the beams' rest frame, put the beams back to back on z, then apply a fixed
  generic rotation.
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

/*
  YFS: VIRTUAL_CHECK's second frame: BeamRestFrameOnZ under a fixed generic
  rotation different from CanonicalBeamFrame's.
*/
inline Vec4D_Vector VirtualCheckFrame(const Vec4D_Vector &p)
{
  Vec4D_Vector q(BeamRestFrameOnZ(p));
  static const double ca(cos(1.9)), sa(sin(1.9)), cb(cos(2.3)), sb(sin(2.3)),
                      cc(cos(0.4)), sc(sin(0.4));
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
