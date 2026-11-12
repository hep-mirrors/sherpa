#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/NLO_Base.H"
#include "YFS/NLO/NLO_Base_Internal.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include <functional>
#include <map>
#include <cstdlib>
#include <iostream>
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include "YFS/NLO/Photon_Counterterm.H"
#include "MODEL/Main/Model_Base.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

/*
  The (n+1)-body point for beta_1(k_j) of final-state photons (REAL_FSR_MAP
  rescale_all). The final legs p[2..] before emission sum to Q, mass sqrt(s').
  Q is kept: in its rest frame the legs keep their directions and are scaled
  by one factor xi to mass sqrt((Q-K)^2), then boosted to Q-K, and the beams
  are rebuilt at sqrt(s'), so the real ME and flux are not taken above s'.
*/
bool NLO_Base::MapMomentaFSR(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || p.size() < 4) return false;
  for (const Vec4D &kj : k) if (!PhotonIsFSR(kj)) return false;
  m_map_reduced = false;
  /*
    Start from the post-emission legs so that the photon-lepton angle (and
    hence the collinear structure of |M_1|^2) is the event's. The other FSR
    photons are absorbed into the pair below.
  */
  Vec4D Kall;
  for (const Vec4D &kj : m_FSRPhotons) Kall += kj;
  /*
    REAL_FSR_MAP pre_emission: the single-emission point, i.e. the
    pre-emission legs with only photon j's recoil applied. It tends to the
    crude configuration as k_j -> 0, so the bracket (r_j - S~ B)/S~ vanishes.
  */
  if (RealFSRMap() == realfsrmap::pre_emission && m_plab.size() == p.size()) {
    for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
    Kall = Vec4D();
  } else if (m_postlab.size() == p.size()) {
    for (size_t i = 2; i < p.size(); ++i) p[i] = m_postlab[i];
  } else Kall = Vec4D();
  Vec4D Q;
  for (size_t i = 2; i < p.size(); ++i) Q += p[i];
  Q += Kall;
  const double sq(Q.Abs2());
  if (!(sq > 0.)) return false;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostQ(Q);
  for (size_t i = 0; i < p.size(); ++i) { pRot.RotateBack(p[i]); boostQ.Boost(p[i]); }
  Vec4D K;
  for (Vec4D &kj : k) { pRot.RotateBack(kj); boostQ.Boost(kj); K += kj; }
  const double M(sqrt(sq));
  const Vec4D Qp(Vec4D(M, 0., 0., 0.) - K);
  const double Mp2(Qp.Abs2());
  if (!(Mp2 > 0.)) return false;
  const double Mp(sqrt(Mp2));
  // Work in the rest frame of the legs' own sum L, where their 3-momenta
  // cancel and a common rescaling keeps them cancelling, then boost onto Q - K.
  Vec4D L;
  for (size_t i = 2; i < p.size(); ++i) L += p[i];
  if (!(L.Abs2() > 0.)) return false;
  Poincare boostL(L), boostQp(Qp);
  for (size_t i = 2; i < p.size(); ++i) boostL.Boost(p[i]);
  // xi: sum_i sqrt(m_i^2 + xi^2 |p_i|^2) = Mp, by bisection (monotone in xi)
  std::vector<double> m2, pp2;
  double msum(0.);
  for (size_t i = 2; i < p.size(); ++i) {
    const double mi(m_flavs[i].Mass());
    m2.push_back(mi*mi); pp2.push_back(Vec3D(p[i]).Sqr()); msum += mi;
  }
  if (msum >= Mp) return false;
  auto etot = [&](double xi) { double e(0.);
    for (size_t j = 0; j < m2.size(); ++j) e += sqrt(m2[j] + xi*xi*pp2[j]);
    return e; };
  double lo(0.), hi(1.);
  while (etot(hi) < Mp) hi *= 2.;
  for (int it = 0; it < 200; ++it) {
    const double mid(0.5*(lo+hi));
    (etot(mid) < Mp ? lo : hi) = mid;
  }
  const double xi(0.5*(lo+hi));
  for (size_t i = 2; i < p.size(); ++i) {
    const Vec3D v(xi*Vec3D(p[i]));
    Vec4D f(sqrt(m2[i-2] + v.Sqr()), v);
    boostQp.BoostBack(f);
    p[i] = f;
  }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double m1 = m_flavs[0].Mass(), m2b = m_flavs[1].Mass();
  const double lamCM = 0.5 * sqrt(Lambda(sq, m1*m1, m2b*m2b) / sq);
  p[0] = {sqrt(lamCM*lamCM + m1*m1), 0, 0,  sign_z * lamCM};
  p[1] = {sqrt(lamCM*lamCM + m2b*m2b), 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) { pRot2.Rotate(p[i]); boostLab.BoostBack(p[i]); }
  for (Vec4D &kj : k) { pRot2.Rotate(kj); boostLab.BoostBack(kj); }
  return true;
}

/*
  The (n+1)-body point for beta_1(k_j) of final-state photons when more than
  one dipole radiates (REAL_FSR_MAP dipole, the default).

  Each FSR photon is radiated by one resonant pair (Photon::Dip()) and only
  that pair recoils, Q_D(pre) = q_1' + q_2' + K_D. Every other pair is put
  back at its pre-emission momentum (m_plab), and k_j's pair is rebuilt at
  Q_D - k_j with its post-emission direction in its own rest frame. Several
  selected photons are grouped by dipole. With one radiating pair this is
  MapMomentaFSR. The beams are rebuilt at sqrt(s') as there.
*/
/*
  The pole scheme's point for photons radiated by the W-pair dipole. The
  dipole names the charged leptons but the W's recoil: the W's are rebuilt
  back to back in the frame of W- + W+ + K and their daughters carried over
  (inverse of Define_Dipoles::ApplyPoleRecoil).
*/
bool NLO_Base::MapMomentaFSRPole(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (p_dipoles == nullptr || !p_dipoles->PoleActive()) return false;
  const YFS::DipoleSet::WWLegs &w(p_dipoles->WW());
  if (!w.ok || w.lm >= p.size() || w.lp >= p.size() || w.nm >= p.size()
      || w.np >= p.size()) return false;
  Vec4D K;
  for (const Vec4D &kj : k) {
    const YFS::Photon *g(FindPhoton(kj));
    if (g == nullptr || !g->IsFSR() || g->Dip() == nullptr) return false;
    if (!g->Dip()->GetFlav(0).IsVector()) return false;   // not the W pair
    K += kj;
  }
  m_map_reduced = false;
  for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
  const std::size_t dau[2][2] = {{w.lm, w.nm}, {w.lp, w.np}};
  Vec4D wpost[2], wpre[2];
  for (int i(0); i < 2; ++i) wpost[i] = m_postlab[dau[i][0]] + m_postlab[dau[i][1]];
  const Vec4D Qp(wpost[0] + wpost[1] + K);
  const double Mp2(Qp.Abs2());
  if (!(Mp2 > 0.) || !(Qp[0] > 0.)) return false;
  const double m1(wpost[0].Mass()), m2(wpost[1].Mass());
  if (!(m1 > 0.) || !(m2 > 0.) || m1 + m2 >= sqrt(Mp2)) return false;
  // the W direction: the post-emission W's in their own rest frame
  Vec4D q1(wpost[0]), q2(wpost[1]);
  const Vec4D L(q1 + q2);
  if (!(L.Abs2() > 0.) || !(L[0] > 0.)) return false;
  Poincare boostL(L);
  boostL.Boost(q1); boostL.Boost(q2);
  Vec3D n(Vec3D(q1) - Vec3D(q2));
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(Mp2, m1*m1, m2*m2)/Mp2));
  wpre[0] = Vec4D(sqrt(m1*m1 + pcm*pcm),  pcm*n);
  wpre[1] = Vec4D(sqrt(m2*m2 + pcm*pcm), -pcm*n);
  Poincare boostQp(Qp);
  boostQp.BoostBack(wpre[0]); boostQp.BoostBack(wpre[1]);
  for (int i(0); i < 2; ++i) {
    Poincare toRest(wpost[i]), fromPre(wpre[i]);
    for (int j(0); j < 2; ++j) {
      Vec4D q(m_postlab[dau[i][j]]);
      toRest.Boost(q);
      fromPre.BoostBack(q);
      p[dau[i][j]] = q;
    }
  }
  return RebuildBeamsAtPreFSR(p, k, "FSRMAPP");
}

bool NLO_Base::MapMomentaFSRDipole(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || p.size() < 4) return false;
  if (m_postlab.size() != p.size() || m_plab.size() != p.size()) return false;
  if (p_dipoles && p_dipoles->PoleActive()) return MapMomentaFSRPole(p, k);
  // The selected photons, grouped by the dipole that radiated them.
  std::map<std::pair<int,int>, Vec4D> ksel;
  for (const Vec4D &kj : k) {
    const YFS::Photon *g(FindPhoton(kj));
    if (g == nullptr || !g->IsFSR() || g->Dip() == nullptr) return false;
    const int l(g->Dip()->Left()), r(g->Dip()->Right());
    if (l < 2 || r < 2 || l >= (int)p.size() || r >= (int)p.size() || l == r)
      return false;
    ksel[std::make_pair(std::min(l,r), std::max(l,r))] += kj;
  }
  m_map_reduced = false;
  // Every final leg at its pre-emission momentum: the pairs that did not
  // radiate a selected photon have their own photons re-absorbed by this.
  for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
  for (const auto &sel : ksel) {
    const int l(sel.first.first), r(sel.first.second);
    const Vec4D K(sel.second);
    const Vec4D Qd(m_plab[l] + m_plab[r]);     // the pair before it radiated
    const Vec4D Qp(Qd - K);                    // the pair with only K taken out
    const double Mp2(Qp.Abs2());
    if (!(Mp2 > 0.) || !(Qp[0] > 0.)) return false;
    const double Mp(sqrt(Mp2));
    const double m1(m_flavs[l].Mass()), m2(m_flavs[r].Mass());
    if (m1 + m2 >= Mp) return false;
    // The decay direction the event has: the post-emission legs in their
    // own rest frame are back to back along it.
    Vec4D q1(m_postlab[l]), q2(m_postlab[r]);
    const Vec4D L(q1 + q2);
    if (!(L.Abs2() > 0.) || !(L[0] > 0.)) return false;
    Poincare boostL(L);
    boostL.Boost(q1); boostL.Boost(q2);
    Vec3D n(Vec3D(q1) - Vec3D(q2));
    if (!(n.Abs() > 0.)) return false;
    n = n/n.Abs();
    const double pcm(0.5*sqrt(Lambda(Mp2, m1*m1, m2*m2)/Mp2));
    Vec4D f1(sqrt(m1*m1 + pcm*pcm),  pcm*n);
    Vec4D f2(sqrt(m2*m2 + pcm*pcm), -pcm*n);
    Poincare boostQp(Qp);
    boostQp.BoostBack(f1); boostQp.BoostBack(f2);
    p[l] = f1; p[r] = f2;
  }
  return RebuildBeamsAtPreFSR(p, k, "FSRMAPD");
}

// The beams at sqrt(Q^2), Q the pre-emission final state (= s'), in the
// same frame convention as MapMomentaFSR: into the Q rest frame, the beams
// along the Born axis there, and back. Shared by the flat and the pole map.
bool NLO_Base::RebuildBeamsAtPreFSR(Vec4D_Vector &p, Vec4D_Vector &k,
                                    const char *tag) {
  Vec4D Q;
  for (size_t i = 2; i < p.size(); ++i) Q += m_plab[i];
  const double sq(Q.Abs2());
  if (!(sq > 0.)) return false;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostQ(Q);
  for (size_t i = 2; i < p.size(); ++i) { pRot.RotateBack(p[i]); boostQ.Boost(p[i]); }
  for (Vec4D &kj : k) { pRot.RotateBack(kj); boostQ.Boost(kj); }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double mb1 = m_flavs[0].Mass(), mb2 = m_flavs[1].Mass();
  const double lamCM = 0.5 * sqrt(Lambda(sq, mb1*mb1, mb2*mb2) / sq);
  p[0] = {sqrt(lamCM*lamCM + mb1*mb1), 0, 0,  sign_z * lamCM};
  p[1] = {sqrt(lamCM*lamCM + mb2*mb2), 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) { pRot2.Rotate(p[i]); boostLab.BoostBack(p[i]); }
  for (Vec4D &kj : k) { pRot2.Rotate(kj); boostLab.BoostBack(kj); }
  return true;
}

/*
  The (n+1)-body point for beta_1(k_j) of an initial-state photon when the
  event has other photons (REAL_MAP beam_axis).

  The weight is beta_1(k_j)/S~(k_j) with S~ the eikonal of the event, so
  |M_1|^2 at X_j needs the same collinear structure. The massless eikonal is
  invariant under rescaling of the legs and beam-axis boosts but not under a
  rotation of k relative to the axis, so the rest-frame map can misplace a
  photon relative to the rebuilt beams when a hard spectator is present.

  Here k_j and the beam directions stay as in the lab, the other photons are
  removed from the beams by their longitudinal projection, one factor lambda
  fixes (lambda P~ - k_j)^2 = s', and the final state is boosted rigidly
  from Q_lab to Q' = P' - k_j. At n = 1 the point is the event itself.
  Final-state photons keep the rest-frame path.
*/
bool NLO_Base::MapMomentaBeamAxis(Vec4D_Vector &p, const Vec4D_Vector &k) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  if (!(sp > 0.)) return false;
  // Everything the event radiated that is not among the selected photons:
  // by conservation, so hidden or unlisted photons are removed as well.
  const Vec4D K(P - Q - ksel);
  const double Kn(Vec3D(K)*n);
  const Vec4D Pt(P - Vec4D(K[0], Kn*n));
  const double Pt2(Pt.Abs2());
  if (!(Pt2 > 0.) || !(Pt[0] > 0.)) return false;
  // lambda^2 Pt^2 - 2 lambda Pt.ksel - (sp - ksel^2) = 0, the positive root.
  const double b(Pt*ksel), c(sp - ksel.Abs2());
  const double disc(b*b + Pt2*c);
  if (!(disc >= 0.)) return false;
  const double lambda((b + sqrt(disc))/Pt2);
  if (!(lambda > 0.) || IsBad(lambda)) return false;
  const Vec4D Pp(lambda*Pt);
  const Vec4D Qp(Pp - ksel);
  if (!(Qp[0] > 0.) || !IsEqual(Qp.Abs2(), sp, 1e-8)) return false;
  // On-shell beams back-to-back along +-n in the P' rest frame; P' has no
  // transverse momentum, so that frame is a longitudinal boost away and the
  // axis survives the boost back.
  const double sj(Pp.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return false;
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toPp(Pp);
  toPp.BoostBack(pa); toPp.BoostBack(pb);
  // The final state, rigidly, from its lab total Q to Q'.
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - ksel);
  const double scale(Max(Pp[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << " in the beam-axis reduction, using the rest-frame map."
                  << std::endl;
      return false;
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  m_map_reduced = (K[0] > 1e-12*P[0]);
  return true;
}

/*
  The KKMC convention for beta_1(k_j) in an n-photon event (REAL_MAP scaled).

  KKMC takes a = k.p_a/(p_a.p_b), b = k.p_b/(p_a.p_b) with the full beams and
  B at the event's s'. This is the n = 1 configuration {P_a, P_b, k_j, Q}
  scaled by x = sqrt(s'/(P - K_sel)^2) <= 1: a and b are scale invariant, the
  pair sits at s', and the final state is boosted rigidly from Q_lab to
  Q' = x(P - k_j). The photon in the point is x k_j, so the eikonal and
  |M_1|^2 both carry 1/x^2: CalculateReal takes the ratio |M_1|^2 flux/(S~ B)
  at the point and multiplies the event's own eikonal (m_map_reduced).
  (Q + k_j)^2 is not preserved, so final-state photons keep the rest-frame
  construction.
*/
bool NLO_Base::MapMomentaScaled(Vec4D_Vector &p, Vec4D_Vector &k, bool anylabel) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  if (!anylabel) for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  const Vec4D R(P - ksel);
  const double R2(R.Abs2());
  if (!(sp > 0.) || !(R2 > 0.) || !(R[0] > 0.)) return false;
  const Vec4D K(P - Q - ksel);
  // x <= 1 up to rounding: (P - K_sel)^2 = (Q + K_others)^2 >= Q^2.
  const double x(Min(1., sqrt(sp/R2)));
  if (!(x > 0.) || IsBad(x)) return false;
  const double sj(x*x*P.Abs2());
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return false;
  // Beams x P_a, x P_b rebuilt on shell along +-n in the P rest frame, then
  // returned to the lab (the identity for balanced beams).
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toP(P);
  toP.BoostBack(pa); toP.BoostBack(pb);
  const Vec4D Qp(pa + pb - x*ksel);
  if (!(Qp[0] > 0.) || !(Qp.Abs2() > 0.)) return false;
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - x*ksel);
  const double scale(Max(P[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << " in the scaled reduction, using the rest-frame map."
                  << std::endl;
      return false;
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  for (Vec4D &kj : k) kj = x*kj;
  m_map_reduced = (K[0] > 1e-12*P[0]);
  return true;
}

/*
  REAL_MAP invariant. Keeps the photon's fractions relative to the full
  beams as in the scaled map and also (Q + k_j)^2, the invariant of the
  final-state emission diagrams. Unknowns: x_a, x_b (directions fixed) and a
  photon scale kappa. Conditions:
      P'^2 = (Q + k_j)^2,   (P' - kappa k_j)^2 = s',
      log kappa = w log x_a + (1-w) log x_b,  w = (1 + cos theta_k)/2.
  Solved by iteration on kappa; of the two roots for P' along the axis the
  one with smaller longitudinal momentum is continuous with P' = P at n = 1.
  With no real solution the scaled construction is used. Not the default
  because it carries a heavy tail for ISR-labelled events near the Z pole.
*/
bool NLO_Base::MapMomentaInvariant(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || m_bornMomenta.size() < 2 || p.size() < 3) return false;
  for (const Vec4D &kj : k) if (PhotonIsFSR(kj)) return false;
  const Vec4D P(m_bornMomenta[0] + m_bornMomenta[1]);
  const Vec3D pa3(m_bornMomenta[0]);
  if (!(pa3.Abs() > 0.)) return false;
  const Vec3D n(pa3/pa3.Abs());
  Vec4D Q, ksel;
  for (size_t i(2); i < p.size(); ++i) Q += p[i];
  for (const Vec4D &kj : k) ksel += kj;
  const double sp(Q.Abs2());
  const Vec4D R(Q + ksel);
  const double sj(R.Abs2());
  if (!(sp > 0.) || !(sj > sp) || !(R[0] > 0.)) return false;
  const Vec4D K(P - R);
  if (!(K[0] > 1e-12*P[0])) return MapMomentaScaled(p, k);   // n = 1: the event
  const double Es(ksel[0]), ks(Vec3D(ksel)*n), D(Es*Es - ks*ks), ksel2(ksel.Abs2());
  if (!(Es > 0.) || !(D > 0.)) return false;
  const double ck(Vec3D(ksel).Abs() > 0. ? ks/Vec3D(ksel).Abs() : 0.);
  const double w(0.5*(1. + ck));
  const double Ea(m_bornMomenta[0][0]), Eb(m_bornMomenta[1][0]);
  const double pza(Vec3D(m_bornMomenta[0])*n), pzb(Vec3D(m_bornMomenta[1])*n);
  // E' = xa Ea + xb Eb, pz' = xa pza + xb pzb (pzb < 0 for balanced beams).
  const double det(Ea*pzb - Eb*pza);
  if (IsZero(det)) return false;
  double kappa(sqrt(sp/(P - ksel).Abs2())), xa(1.), xb(1.), Ep(0.), pzp(0.);
  bool ok(false);
  for (int it(0); it < 50; ++it) {
    const double C((sj + kappa*kappa*ksel2 - sp)/(2.*kappa));   // P'.ksel
    const double disc(C*C - D*sj);
    if (!(disc >= 0.)) break;
    const double r1((C*ks + Es*sqrt(disc))/D), r2((C*ks - Es*sqrt(disc))/D);
    pzp = dabs(r1) < dabs(r2) ? r1 : r2;
    Ep  = (C + pzp*ks)/Es;
    if (!(Ep > 0.) || !(Ep*Ep - pzp*pzp > 0.)) break;
    // P' = xa P_a + xb P_b along the axis.
    xa = (Ep*pzb - Eb*pzp)/det;
    xb = (Ea*pzp - pza*Ep)/det;
    if (!(xa > 0.) || !(xb > 0.)) break;
    const double kn(exp(w*log(xa) + (1. - w)*log(xb)));
    if (dabs(kn - kappa) < 1e-12*kn) { kappa = kn; ok = true; break; }
    kappa = kn;
  }
  if (!ok || IsBad(kappa) || !(kappa > 0.)) return MapMomentaScaled(p, k);
  const Vec4D Pp(Ep, pzp*n);
  const double m1(m_flavs[0].Mass()), m2(m_flavs[1].Mass());
  if (sj <= sqr(m1 + m2)) return MapMomentaScaled(p, k);
  const double lamCM(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(lamCM*lamCM + m1*m1),  lamCM*n);
  Vec4D pb(sqrt(lamCM*lamCM + m2*m2), -lamCM*n);
  Poincare toPp(Pp);
  toPp.BoostBack(pa); toPp.BoostBack(pb);
  const Vec4D Qp(pa + pb - kappa*ksel);
  if (!(Qp[0] > 0.) || !IsEqual(Qp.Abs2(), sp, 1e-6)) return MapMomentaScaled(p, k);
  Poincare fromQ(Q), toQp(Qp);
  Vec4D_Vector q(p.begin() + 2, p.end());
  Vec4D qsum;
  for (Vec4D &qi : q) { fromQ.Boost(qi); toQp.BoostBack(qi); qsum += qi; }
  Vec4D bal(pa + pb - qsum - kappa*ksel);
  const double scale(Max(P[0], 1.));
  for (int mu(0); mu < 4; ++mu)
    if (dabs(bal[mu]) > 1e-7*scale) {
      msg_Error() << METHOD << "(): momentum imbalance " << bal
                  << ", using the scaled reduction." << std::endl;
      return MapMomentaScaled(p, k);
    }
  p[0] = pa; p[1] = pb;
  for (size_t i(2); i < p.size(); ++i) p[i] = q[i-2];
  for (Vec4D &kj : k) kj = kappa*kj;
  m_map_reduced = true;
  return true;
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D_Vector &k) {
  // YFS: REAL_MAP and REAL_FSR_MAP, documented in YFS_Base::RegisterDefaults
  static const realmap::code mapmode(ATOOLS::Settings::GetMainSettings()["YFS"]
                                     ["REAL_MAP"].Get<realmap::code>());
  m_map_reduced = false;
  if (mapmode == realmap::beam_axis && MapMomentaBeamAxis(p, k)) return;
  if (mapmode == realmap::scaled && MapMomentaScaled(p, k)) return;
  if (mapmode == realmap::invariant && MapMomentaInvariant(p, k)) return;
  const realfsrmap::code fsrmap(RealFSRMap());
  if (fsrmap == realfsrmap::dipole && MapMomentaFSRDipole(p, k)) return;
  // dipole falls back to the common rescaling when the per-dipole
  // construction does not apply (a photon without a dipole, a pair below
  // threshold).
  if (fsrmap != realfsrmap::rest_frame && MapMomentaFSR(p, k)) return;
  if (mapmode != realmap::rest_frame) {
    // Only initial-state photons are reduced by the new constructions; a
    // final-state photon (or an early exit) comes here by design. Counted so
    // that a solver that never converges cannot masquerade as the legacy map.
    static long nfall(0);
    if (++nfall == 1000 && msg_LevelIsDebugging())
      msg_Debugging() << METHOD << "(): 1000 photons on the rest-frame map with REAL_MAP "
                      << mapmode << std::endl;
  }
  Vec4D Q;
  Vec4D QQ;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  for (size_t i = 2; i < p.size(); ++i) Q += p[i];
  Vec4D ksum;
  for (const Vec4D &ki : k) ksum += ki;
  Q += ksum;
  const double sq = Q.Abs2();
  Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  for (size_t i = 0; i < p.size(); ++i) {
    pRot.RotateBack(p[i]);
    boostQ.Boost(p[i]);
  }
  for (Vec4D &ki : k) {
    pRot.RotateBack(ki);
    boostQ.Boost(ki);
  }
  // ksum must be recomputed: the photons have been boosted since.
  ksum = Vec4D();
  for (const Vec4D &ki : k) ksum += ki;
  for (size_t i = 2; i < p.size(); ++i) QQ += p[i];
  QQ += ksum;
  const double sqq = QQ.Abs2();
  if (!IsEqual(sqq, sq, 1e-6)) {
    msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD
                << std::endl;
  }
  const double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  const double m1 = m_flavs[0].Mass();
  const double m2 = m_flavs[1].Mass();
  const double lamRaw = sqq*sqq + sqr(m1*m1) + sqr(m2*m2)
                        - 2.*sqq*m1*m1 - 2.*sqq*m2*m2 - 2.*m1*m1*m2*m2;
  if (lamRaw < 0.)
    msg_Error()<<METHOD<<"(): below-threshold Kaellen argument = "<<lamRaw
               <<" (sqq = "<<sqq<<", threshold = "<<sqr(m1+m2)<<")"<<std::endl;
  const double lamCM = 0.5 * sqrt(Lambda(sqq, m1 * m1, m2 * m2) / sqq);
  const double E1 = lamCM * sqrt(1 + m1 * m1 / sqr(lamCM));
  const double E2 = lamCM * sqrt(1 + m2 * m2 / sqr(lamCM));
  p[0] = {E1, 0, 0, sign_z * lamCM};
  p[1] = {E2, 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (size_t i = 0; i < p.size(); ++i) {
    pRot2.Rotate(p[i]);
    boostLab.BoostBack(p[i]);
  }
  for (Vec4D &ki : k) {
    pRot2.Rotate(ki);
    boostLab.BoostBack(ki);
  }
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
  Vec4D_Vector ks{k};
  MapMomenta(p, ks);
  k = ks[0];
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k,
                          const Vec4D_Vector &ref) {
  const Vec4D_Vector save(m_bornMomenta);
  m_bornMomenta = ref;
  Vec4D_Vector ks{k};
  MapMomenta(p, ks);
  k = ks[0];
  m_bornMomenta = save;
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k1, Vec4D &k2) {
  Vec4D_Vector ks{k1, k2};
  MapMomenta(p, ks);
  k1 = ks[0];
  k2 = ks[1];
}
