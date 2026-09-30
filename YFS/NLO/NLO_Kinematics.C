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

void NLO_Base::RandomRotate(Vec4D &p) {
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta) * t1[2] - sin(m_ranTheta) * t1[3];
  p[3] = sin(m_ranTheta) * t1[2] + cos(m_ranTheta) * t1[3];
  t1 = p;
  // rotate around z
  p[1] = cos(m_ranPhi) * t1[1] - sin(m_ranPhi) * t1[2];
  p[2] = sin(m_ranPhi) * t1[1] + cos(m_ranPhi) * t1[2];
}

void NLO_Base::CheckMappingRecoil(const Vec4D_Vector &p, const Vec4D &ksum) {
  Vec4D q;
  for (size_t i(2); i < p.size(); ++i) q += p[i];
  const Vec3D imbalance(Vec3D(ksum) + Vec3D(q));
  const double scale(Max(Vec3D(ksum).Abs(), Vec3D(q).Abs()));
  if (scale <= 0.) return;
  const double rel(imbalance.Abs() / scale);
  if (rel > 1e-5)
    msg_Error() << METHOD << "(): YFS mapping recoil imbalance " << rel
                << " (photons " << Vec3D(ksum) << ", hard system " << Vec3D(q)
                << ")" << std::endl;
}

/*
  The (n+1)-body point for beta_1(k_j) of FINAL-state photons (REAL_FSR_MAP: rescale_all).

  p[2..] are the final legs BEFORE final-state emission (m_reallab), so their
  sum is the ISR-reduced total momentum Q, of mass sqrt(s'). The rest-frame
  map below adds k on top and rebuilds the beams at sqrt((Q+k)^2), i.e. ABOVE
  sqrt(s'): e+e- -> mu mu at 0.7 GeV gave sqrt(s_j) = 0.83-0.89 GeV for hard
  wide-angle FSR photons, the real ME and flux were taken at that unphysical
  point while the subtraction used S~ B at s', and the fixed-order photon
  spectrum grew a heavy-weight bump above E_gamma ~ 0.12 GeV (independent of
  the ME provider and of REAL_MAP, which touches ISR photons only).

  Here Q is kept: in the Q rest frame the final legs keep their directions and
  are rescaled by one factor xi so that they carry the mass sqrt((Q-K)^2),
  then boosted to move with Q-K; the beams are rebuilt at sqrt(Q^2) = sqrt(s').
  With a single FSR photon this is the YFS final-state recoil of the event.
*/
bool NLO_Base::MapMomentaFSR(Vec4D_Vector &p, Vec4D_Vector &k) {
  if (k.empty() || p.size() < 4) return false;
  for (const Vec4D &kj : k) if (!PhotonIsFSR(kj)) return false;
  m_map_reduced = false;
  /*
    Start from the legs AFTER final-state emission, not before. With the
    pre-emission legs the photon-lepton angle of the point is not the event's,
    and the collinear structure of |M_1|^2 is off by large factors: on e+e- ->
    mu mu at 0.7 GeV, single-FSR-photon events had r/rho_1(CEEX) spread over
    0.003-0.5 where single-ISR-photon events give 0.01613 on every event. The
    other FSR photons are absorbed into the pair below (their sum is added to
    the target mass); with one FSR photon the point IS the event, as for CEEX.
  */
  Vec4D Kall;
  for (const Vec4D &kj : m_FSRPhotons) Kall += kj;
  /*
    REAL_FSR_MAP: pre_emission (3) - the SINGLE-EMISSION point: the pre-emission legs with
    photon j's own recoil applied by the rescaling below and nothing else.
    The post-emission construction above re-absorbs the other photons by
    rescaling, which keeps the post-emission DIRECTIONS; a soft photon whose
    parent lepton was kicked by a hard companion then sits at a point whose
    collinear structure is that of the post-emission legs while the crude it
    is subtracted against (REAL_SUB_EIK 3, the generation density) is the
    pre-emission one, and its bracket (r_j - S~ B)/S~ is O(1) instead of
    vanishing: Z-pole mu mu, nfsr >= 3 all soft, 10th percentile of the
    fixed-order factor -5. Starting from the pre-emission legs the point
    tends to the crude configuration as k_j -> 0, and with one photon it is
    the generator's own post-emission point (same rescaling recipe).
  */
  { if (RealFSRMap() == realfsrmap::pre_emission && m_plab.size() == p.size()) {
      for (size_t i = 2; i < p.size(); ++i) p[i] = m_plab[i];
      Kall = Vec4D();
    } else if (m_postlab.size() == p.size()) {
      for (size_t i = 2; i < p.size(); ++i) p[i] = m_postlab[i];
    } else Kall = Vec4D(); }
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
  // Into the rest frame of the legs' own sum L, where their 3-momenta cancel
  // and a common rescaling keeps them cancelling; the rescaled legs are then
  // boosted onto Q - K. (Boosting into the Q - K frame instead left the
  // absorbed photons' 3-momentum unbalanced whenever K != all FSR photons.)
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
  { Vec4D res(p[0] + p[1]);
    for (size_t i = 2; i < p.size(); ++i) res -= p[i];
    for (const Vec4D &kj : k) res -= kj;
    if (Vec3D(res).Abs() > 1e-9 || std::abs(res[0]) > 1e-9) {
      Vec4D post; for (size_t i = 2; i < m_postlab.size(); ++i) post += m_postlab[i];
      Vec4D pre;  for (size_t i = 2; i < m_plab.size(); ++i)    pre  += m_plab[i];
      Vec4D isr;  for (const Vec4D &kj : m_ISRPhotons) isr += kj;
      std::cerr<<std::setprecision(10)<<"@@@ FSRMAP residual="<<res
               <<" xi="<<xi<<" Mp="<<Mp<<" msum="<<msum
               <<" nFSR="<<m_FSRPhotons.size()<<" nk="<<k.size()
               <<"\n   Kall(lab)="<<Kall<<" post(lab)="<<post<<" pre(lab)="<<pre
               <<" isr(lab)="<<isr<<" P="<<(m_bornMomenta[0]+m_bornMomenta[1])
               <<"\n   post+Kall+isr-P="<<(post+Kall+isr-m_bornMomenta[0]-m_bornMomenta[1])
               <<" pre+isr-P="<<(pre+isr-m_bornMomenta[0]-m_bornMomenta[1])<<std::endl;
    } }
  return true;
}

/*
  The (n+1)-body point for beta_1(k_j) of FINAL-state photons when the final
  state has MORE THAN ONE radiating dipole (REAL_FSR_MAP: dipole, the default).

  Each final-state photon is radiated by one resonant pair (its dipole,
  YFS::Photon::Dip()): Dipole::GenerateEmissions samples it in that pair's
  rest frame and YFS_Handler::CalculateFSR lets only that pair recoil, so in
  the event every pair separately satisfies  Q_D(pre) = q_1' + q_2' + K_D.
  MapMomentaFSR above absorbs the OTHER photons by rescaling ALL final legs
  together: with two pairs the pair that did not radiate k_j is moved off its
  pre-emission momentum and k_j's pair does not land on Q_D - k_j, so both
  resonance propagators of |M_1|^2 are evaluated off the event's invariants
  while the subtraction S~ B sits on the event's ones. On the Z pole a shift
  of a fraction of Gamma_Z is a large factor, which is what the numbers show.
  e+e- -> mu mu tau tau at 250 GeV (doc/examples/YFS/zpole/br_table/
  250_eemumu, 8k events, mu-pair mass within 2 GeV of M_Z):
      nfsr = 1: <BR>/<CEEX> = 0.987, median BR/CEEX = 1.06   (the event)
      nfsr = 2: <BR>/<CEEX> = 0.21,  median 1.03 - the mean is a tail;
      nfsr = 3: 0.17;   nfsr >= 4: 0.12;
      a hard FSR photon (x > 0.05) with other FSR photons: median BR = -0.66,
      <BR> = -3.5 against <CEEX> = 0.75.
  In the Rivet Z1_mass overlay the Born+real dipped 25% at the mu-pair peak
  where LO, EEX and CEEX agree in shape.

  Here every pair other than k_j's is put back at its pre-emission momentum
  (m_plab, its photons re-absorbed), and k_j's pair is rebuilt at
  Q_D(pre) - k_j: its two legs keep the direction they have in the event in
  their own rest frame (the post-emission legs of m_postlab boosted to rest,
  which is the direction the generator drew), are given the momentum of a
  two-body decay of mass sqrt((Q_D - k_j)^2), and are boosted onto Q_D - k_j.
  Several selected photons (CalculateRealReal) are grouped by dipole and
  absorbed by their own pair each. Legs in no radiating dipole are untouched.
  When k_j is the only photon of its pair the pair is the event's, and with
  one radiating pair in the process the whole construction is MapMomentaFSR.
  The beams are rebuilt at sqrt(Q^2) = sqrt(s') exactly as there.
*/
/*
  The pole scheme's (n+1)-body point for photons the W-PAIR dipole radiated.
  Such a photon's dipole names the two charged leptons (Left/Right are the
  lepton positions, DipoleSet::BuildPole), but the pair that recoiled is the
  W's, so the flat map above this would rebuild the lepton pair, the wrong
  pre-emission point by the W recoil. Here: the W pair before the selected
  photons, Qp = W- + W+ + K, the two W's back to back in Qp's frame along
  their post-emission direction with their (preserved) masses, and each W's
  daughters carried from the post-emission W to the rebuilt one - the
  inverse of Define_Dipoles::ApplyPoleRecoil. Photons from any other dipole
  (none exist in the pole scheme) fall back to the flat map.
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
  { Vec4D res(p[0] + p[1]);
    for (size_t i = 2; i < p.size(); ++i) res -= p[i];
    for (const Vec4D &kj : k) res -= kj;
    const double scale(Max(1., p[0][0]));
    if (Vec3D(res).Abs() > 1e-9*scale || std::abs(res[0]) > 1e-9*scale) {
      static long nprint(0);
      if (++nprint <= 20)
        std::cerr<<std::setprecision(10)<<"@@@ "<<tag<<" residual="<<res
                 <<" nk="<<k.size()
                 <<" nFSR="<<m_FSRPhotons.size()<<" P="<<(m_bornMomenta[0]+m_bornMomenta[1])
                 <<" Q="<<Q<<std::endl;
    } }
  return true;
}

/*
  The (n+1)-body point at which beta_1(k_j) is evaluated when the event has
  OTHER photons, for a photon radiated from the initial state (REAL_MAP: 1).

  What the point has to reproduce is the soft-photon radiation pattern of
  the event: the weight is beta_1(k_j)/S~(k_j) with S~ the eikonal of the
  EVENT, so |M_1|^2 at X_j must carry the same collinear structure for k_j
  or the ratio is unbounded. The massless eikonal sum_i Q_i p_i/(p_i.k) is
  invariant under a rescaling p_i -> x_i p_i of every leg and, for the beams,
  under boosts along the beam axis (it is 4/k_perp^2), but NOT under a
  rotation or a transverse boost of k relative to the axis.

  The rest-frame construction below (REAL_MAP: 0) boosts final state + k_j
  into their common rest frame and rebuilds the beams along the lab axis
  THERE. With a hard spectator photon that frame moves at beta ~ 0.7 and the
  aberrated k_j can land on the rebuilt beam axis: measured on a 250 GeV
  radiative-return event, a 7.4 GeV photon at 34 degrees to the beams (lab
  eikonal 4.3e-5) came out at 20 mrad (eikonal 1.1e-2), and a -21% real
  correction became -52 times the Born. The tail of that ratio is 1/f, so
  the cross-section error never shrank with statistics (2.5-3.5% at 100k for
  250 GeV nu nu, mu mu, tau tau; per-event factors down to -2600).

  Here instead:
    - k_j and the beam DIRECTIONS stay as in the lab;
    - the other photons are taken out of the beams through their
      longitudinal projection, P~ = P - (E_K, (K.n) n), which is the
      light-cone removal that leaves x_a P_a + x_b P_b on the axis;
    - one common factor lambda fixes (lambda P~ - k_j)^2 = s', so the pair
      invariant - the Z propagator - is exactly the event's and matches the
      Born of the subtraction;
    - the final state is boosted rigidly from Q_lab to Q' = P' - k_j, so its
      internal invariants and masses are untouched.
  The II eikonal is then preserved up to the electron-mass terms (whose
  ratio lies in [x^4, 1]: the dead cone only shrinks the term), and at n = 1
  there is nothing to remove, lambda = 1, and the point is the event itself,
  so the one-photon stream is unchanged. What is given up is the angle of
  k_j to the final-state fermions, which changes by the bounded Doppler
  factor of the Q_lab -> Q' boost; for a photon radiated from the FINAL
  state that is the singular structure, so those keep the rest-frame path.
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
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" lambda="<<lambda<<" beamaxis=1"<<std::endl; }
  return true;
}

/*
  The KKMC convention for beta_1(k_j) in an n-photon event (REAL_MAP: 2).

  KKMC's EEX beta_1 for an initial-state photon is S~(k) B(s') times
  [(1-a)^2 + (1-b)^2]/2 - 1 with a = k.p_a/(p_a.p_b), b = k.p_b/(p_a.p_b)
  taken with the FULL beams and B at the event's s': the other photons enter
  only through s'. The beam-axis reduction above judges the photon against
  beams already reduced by the others, a' = a/x_b, which for two hard photons
  counts the energy loss twice. Leading-log check, two collinear ISR photons
  with fractions x_1, x_2 of the beam, exact factor
  (1+(1-x_1)^2)/2 * (1+(1-x_2/(1-x_1))^2)/2 symmetrised:
      x_1 = x_2 = 0.44 : exact 0.34, KKMC 0.32, reduced beams 0.04
      x_1 = 0.8, x_2 = 0.05 : exact 0.45, KKMC 0.47, reduced beams 0.29
  so the KKMC convention is the better O(alpha) truncation (the remainder is
  the genuine beta_2), and it is what YFS.EEX and CEEX effectively carry.
  Measured, 250 GeV nu nu, 50k events, same sample: the beam-axis reduction
  gives YFS.BR 2.838 +- 0.3%, this construction 3.166 +- 0.2%, CEEX (no
  virtual) 3.272 +- 0.4%, EEX 3.513 (with its +7% virtual); the rest-frame
  map gave 2.925 +- 10%. The remaining 3% to CEEX sits in the events where
  two hard photons share the radiative return, where CEEX's coherent
  amplitude-level sum carries an effective beta_2 that O(alpha) cannot.

  A momentum-conserving point with the full-beam a, b and the pair at s' is
  the n = 1 configuration {P_a, P_b, k_j, Q} scaled by one factor x:
      x^2 (P - k_j)^2 = s'   ->   x = sqrt(s'/(P - K_sel)^2) <= 1,
  a and b are scale invariant, the pair sits at s' by construction, and the
  final state is boosted rigidly from Q_lab to Q' = x(P - k_j). The photon in
  the point is x k_j, so the eikonal of the point is the event's over x^2 and
  the SAME 1/x^2 sits in |M_1|^2: CalculateReal therefore takes the ratio
  |M_1|^2 flux/(S~ B) at the point and multiplies the event's own eikonal
  (m_map_reduced). At n = 1, x = 1 and nothing changes. What is not preserved
  is (Q + k_j)^2, the invariant of the final-state emission diagrams of the
  full real ME; for photons radiated from the final state that IS the
  resonance, so those keep the rest-frame construction.
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
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" x="<<x<<" scaled=1"<<std::endl; }
  return true;
}

/*
  REAL_MAP: 3. The scaled construction keeps the photon's fractions a, b
  relative to the full beams but moves (Q + k_j)^2, the invariant mass of
  final state plus photon. For a charged final state the full real ME also
  has the diagrams with the photon radiated from the final state, whose Z
  propagator sits at exactly that invariant: with the Z-pole Born, moving it
  by a few GeV changes |M_1|^2 by orders of magnitude. Measured: 250 GeV
  mu mu with REAL_MAP: 2 kept the central value (BR 3.21, CEEX 3.22) but a
  5% error from 466 events with |BR| > 10, and the Z-pole mu mu Born+real
  moved by 0.5%.

  Here both are kept - and the measurement is the reason this is NOT the
  default: 250 GeV mu mu, 100k, REAL_MAP: 3 gives BR 4.98 +- 26% against
  3.21 +- 4.8% and 2.82 +- 5.4% (two versions of the denominator, same
  seeds) for REAL_MAP: 2 (CEEX 3.22) - the 5% is a heavy tail, not a
  Gaussian error, and the central value moves with it. Keeping the invariant exactly
  puts more ISR-labelled events onto the Z pole of their final-state
  emission diagrams, where the full |M_1|^2 is shared between the I and F
  labels in the ratio S~_II/S~_tot but divided by the ISR label's own
  S~_tot B((mu mu)^2), which is orders of magnitude below the physics. The
  label sum is exact (both labels together give the exact real once), so
  this is variance, not bias; the cure is a multi-channel denominator
  S~_II B_I + S~_FF B_F, which needs the Born at the other label's
  kinematics and is not attempted here. At the Z pole modes 2 and 3 agree
  per event to four digits. Unknowns: the beam energies x_a, x_b (directions fixed)
  and one scale kappa for the photon. Conditions:
      P'^2 = (Q + k_j)^2            (final-state-emission invariant),
      (P' - kappa k_j)^2 = s'       (the pair at the event's s'),
      log kappa = w log x_a + (1-w) log x_b,  w = (1 + cos theta_k)/2,
  the last being the scaled convention for the beam the photon is collinear
  with (kappa = x_a exactly along beam a: its fraction of that beam is the
  event's) and a smooth interpolation in between. The photon's energy in the
  (Q + k_j) rest frame is then the event's too, so the final-state emission
  is as hard as it really was; only the angle between photon and fermions
  changes, by the bounded Doppler factor of the Q_lab -> Q' boost. Solved by
  iteration on kappa; the 2x2 system for P' along the axis is quadratic,
  the root with the smaller longitudinal momentum is the one continuous
  with P' = P at n = 1. When it has no real solution (a photon with more
  transverse momentum than energy in the (Q + k_j) frame) the scaled
  construction is used instead.
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
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<sqrt(sj)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksel[0]<<" kappa="<<kappa<<" xa="<<xa<<" xb="<<xb
                     <<" invariant=1"<<std::endl; }
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
  { const realfsrmap::code fsrmap(RealFSRMap());
    if (fsrmap == realfsrmap::dipole && MapMomentaFSRDipole(p, k)) return;
    // dipole falls back to the common rescaling when the per-dipole
    // construction does not apply (a photon without a dipole, a pair below
    // threshold).
    if (fsrmap != realfsrmap::rest_frame && MapMomentaFSR(p, k)) return; }
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
  // CheckMappingRecoil(p, ksum);
  for (size_t i = 2; i < p.size(); ++i) QQ += p[i];
  QQ += ksum;
  const double sqq = QQ.Abs2();
  if (!IsEqual(sqq, sq, 1e-6)) {
    msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD
                << std::endl;
  }
  { static const bool dm(ATOOLS::Settings::GetMainSettings()["YFS"]["REAL_STAB"].Get<int>()!=0);
    if (dm) std::cerr<<"@@@ MAPQ sqrt_sqq="<<(sqq>0?sqrt(sqq):-1.)
                     <<" sqrt_s="<<sqrt(m_s)
                     <<" sqrt_sborn="<<(m_bornMomenta[2]+m_bornMomenta[3]).Mass()
                     <<" Eksum="<<ksum[0]<<std::endl; }
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
