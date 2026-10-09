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
  YFS: RV_MODE - how the squared-level real-virtual (NLO_Part with E) is built
  (name or old integer).

  legacy (0, default): the legacy RV, CalculateRealVirtual(Vec4D): per photon
    [r_RV flux - (B~ rtree flux + S~ B v)/kappa]/S~_crude, the whole
    beta_1^(2) - beta_1^(1) = R v_{n+1} - S~ B v, with its own point, crude,
    flux and coupling, and with E present REAL_COMBINE's product and
    VIRTUAL_COMBINE switched off. IR unsafe: the (n+1) loop runs at
    mu = 1 GeV against the Born's mu^2 = s, and the full-EW finite part is
    mu dependent beyond its IR pole (+0.4% G_mu, +0.8% alpha(0) per event
    whose photons are all soft).

  remainder (1): with v = V_sub/B of the event, v_{n+1} the same object at
    photon j's (n+1)-body point and rho_j = R_j/(crude_j B) the photon's
    exact-over-crude ratio, the exact one-photon O(alpha^2) weight
      beta_0 (1 + v) + [R (1 + v_{n+1}) - S B (1 + v)]/crude
    is identically
      (1 + v)(1 + delta) + rho (v_{n+1} - v),
    delta = (R - S B)/(crude B) the real's bracket. Here v_{n+1} - v is
    taken as lt - V/B, the (n+1) loop over tree minus the Born loop over
    Born, i.e. both with the YFS subtraction on the generation legs.
  remainder_own_legs (2): the same remainder with v_{n+1} the (n+1) loop
    with the YFS subtraction on the (n+1)-body point's OWN legs (post-
    emission legs for a final-state photon), v_{n+1} = lt - B(point legs),
    plus the exponent's leg difference Y_FF+IF(point legs) - Y_FF+IF(
    generation legs) at the exponent's own cutoffs (FormFactorSumFinalLegs).
    That is the subtraction the exponent radiating from the event's legs
    requires: the FSR generator's Sudakov (FSR::YFS_FORM, the log(1/Delta)
    and B-tilde pieces) is on the post-emission legs, and only the handler's
    Kmax-level Y_FF(sqrt(s')/2) and the IF form factor are on the generation
    legs, so v_{n+1} must carry B(point legs) and the Kmax-level difference
    (NOTES-yfsnlo-rr-cutdep-2026-09-30.md sec. 6). The two modes differ by
    2 alpha [Re B(point legs) - Re B(generation legs)]: nothing for an
    initial-state photon (same legs), about 8% of the real-virtual for a
    hard final-state one (Z-pole mu mu, photon tagged). The ACRAIC averaged
    real-virtual (ACRAIC: REAL_VIRTUAL averaged) has always used the own-legs
    v_{n+1} (NLO_Base::m_rvdv).
*/
/*
  YFS: RR_MODE - how the squared-level double real (NLO_Part with W) is built
  (name or old integer).
  legacy (0, default): the legacy assembly (CalculateRealReal(k1,k2), steered by
    RR_CONVENTIONS), which also switches REAL_COMBINE's product off.
  exact (1): exact at O(alpha^2) for every photon pair. f_k is photon k's factor in
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
    NOT for use with IFI_Real false (the paper scheme): there the single
    real's bracket subtracts the coherent eikonal, i.e. the initial-final
    interference of its photon, because the IF form factor at sqrt(s)/2
    already holds it; f_i f_j therefore lacks the IF interference of both
    photons while R_2 has it, and rho_ij - f_i f_j puts the exponentiated IF
    interference back for every pair (O(1) per pair for a wide-angle soft
    photon: Z-pole mu mu, soft pairs, rms 0.4-0.6 B per pair, all noise).
  singles_product (2): the pair's exact over crude beyond the product of the
    singles' exact over crude,
        RR_ij = B (rho_ij - rho_i rho_j) prod_{k != i,j} f_k,
    rho_k = kappa r_k F_k/(c_k B) the single real's own ratio (RRPhotonInfo::rho,
    without its IF subtraction a_k = rho_k - f_k). With two photons the real
    weight is then rho_ij - a_j rho_i - a_i rho_j + a_i a_j: R_2 with the IF
    eikonal of each photon times the other photon's single real taken out,
    which is what the exponentiated IF interference (IFI_Real false) plus the
    one-photon events' brackets already supply; with IFI_Real true a_k = 0
    and this is rho_ij, so the value is right in both schemes. As one photon
    goes soft the remainder tends to rho_i times the difference of the
    coherence ratios S~_coh/S~_crude on the event's legs and on the single's
    point legs, which is O(alpha) times the hard photon's recoil, so the
    single-soft limits hold to the accuracy of the single real's point
    convention; the RR is taken over every generated pair (RR_SOFT_CUT 0),
    so no cut separates a product region from the exact one.
    v x RR and v x R_n are O(alpha^3) and are left out of the squared-level
    weight (YFS_Handler::CalculateNLO); the O(alpha) columns of a W run are
    then those of the same run without W.
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

// YFS: REAL_STAB (diagnostic): dump the real's mapped invariants (DumpMapQ)
inline bool RealStabDump() {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>() != 0);
  return on;
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
    no_flux (1, default since 2026-09-26): none;
    flux_squared (2): the event flux squared, (m_ff^2/s')^2;
    on_crude (3): no flux on r; F = (q + k)^2/q^2 on the crude's final-state
      part instead, q the post-emission pair;
    own_dipole (4): the flux of the photon's own pair, (Q_D - k)^2/Q_D^2.
  The double real takes the pair's ISR flux with no_flux and Flux_Mode event.
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

/*
  YFS: VIRTUAL_CHECK's second frame: BeamRestFrameOnZ under a fixed generic
  rotation different from CanonicalBeamFrame's (z-y-z Euler 1.9, 2.3, 0.4).
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
