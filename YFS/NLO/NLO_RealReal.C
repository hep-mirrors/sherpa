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
#include "YFS/NLO/NLO_Base_Internal.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

/*
  RR_MODE 1: the point pt (legs, then photons) with photon `which` (index
  into pt) taken out and its momentum absorbed by the final state: the final
  legs, in their rest frame, rescaled (masses kept) to the invariant mass
  (Q_f + k)^2 and boosted onto Q_f + k. Beams and the other photons are
  untouched. A soft reduction: as k -> 0 the result tends to pt without k.
*/
static bool AbsorbPhoton(const Vec4D_Vector &pt, size_t nlegs, size_t which,
                         Vec4D_Vector &out)
{
  if (which < nlegs || which >= pt.size() || nlegs < 3) return false;
  const Vec4D k(pt[which]);
  Vec4D Qf;
  for (size_t i(2); i < nlegs; ++i) Qf += pt[i];
  const Vec4D T(Qf + k);
  const double MT2(T.Abs2());
  if (!(MT2 > 0.) || !(Qf.Abs2() > 0.)) return false;
  const double MT(sqrt(MT2));
  Poincare rest(Qf);
  Vec4D_Vector f;
  std::vector<double> m2;
  for (size_t i(2); i < nlegs; ++i) {
    Vec4D v(pt[i]);
    rest.Boost(v);
    f.push_back(v);
    m2.push_back(Max(0., v.Abs2()));
  }
  double msum(0.);
  for (double m : m2) msum += sqrt(m);
  if (!(MT > msum)) return false;
  double xi(1.);
  for (int it(0); it < 60; ++it) {
    double E(0.), dE(0.);
    for (size_t i(0); i < f.size(); ++i) {
      const double p2(Vec3D(f[i]).Sqr()), e(sqrt(xi*xi*p2 + m2[i]));
      E += e;
      if (e > 0.) dE += xi*p2/e;
    }
    if (!(dE > 0.)) return false;
    const double step((E - MT)/dE);
    xi -= step;
    if (std::abs(step) < 1e-15*Max(1., xi)) break;
  }
  if (!(xi > 0.)) return false;
  out.assign(pt.begin(), pt.begin() + 2);
  Poincare onto(T);
  for (size_t i(0); i < f.size(); ++i) {
    const Vec3D p3(xi*Vec3D(f[i]));
    Vec4D v(sqrt(p3.Sqr() + m2[i]), p3);
    onto.BoostBack(v);
    out.push_back(v);
  }
  for (size_t i(nlegs); i < pt.size(); ++i) if (i != which) out.push_back(pt[i]);
  return true;
}

/*
  RR_MODE 1: the point pt with final-state photon `which` recombined with the
  leg of its radiating pair it is closer to in angle (the emitter; the other
  leg is the spectator): in the rest frame of emitter + spectator + photon
  the two legs are rebuilt back to back along emitter + photon, on shell.
  Only the pair changes, and the emitter keeps the direction of
  emitter + photon, so the collinear structure of the OTHER photon is kept
  (absorbing into the whole final state moved the legs by the hard photon's
  recoil, and the eikonals of a collinear companion by factors of 50).
*/
static bool RecombineWithEmitter(const Vec4D_Vector &pt, size_t which, int l, int r,
                                 const ATOOLS::Flavour_Vector &fl, Vec4D_Vector &out)
{
  if (which >= pt.size() || l < 2 || r < 2 || l == r
      || l >= (int)fl.size() || r >= (int)fl.size()) return false;
  const Vec4D k(pt[which]);
  auto angle = [](const Vec4D &a, const Vec4D &b) {
    return Vec3D(a)*Vec3D(b)/(Vec3D(a).Abs()*Vec3D(b).Abs()); };
  const int e(angle(pt[l], k) >= angle(pt[r], k) ? l : r), sp(e == l ? r : l);
  const Vec4D Q(pt[e] + pt[sp] + k);
  const double Q2(Q.Abs2()), me(fl[e].Mass()), ms(fl[sp].Mass());
  if (!(Q2 > sqr(me + ms)) || !(Q[0] > 0.)) return false;
  Poincare rest(Q);
  Vec4D q(pt[e] + k);
  rest.Boost(q);
  Vec3D n(q);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(Q2, me*me, ms*ms)/Q2));
  Vec4D pe(sqrt(me*me + pcm*pcm), pcm*n), ps(sqrt(ms*ms + pcm*pcm), -pcm*n);
  rest.BoostBack(pe);
  rest.BoostBack(ps);
  out.clear();
  for (size_t i(0); i < pt.size(); ++i) {
    if (i == which) continue;
    if ((int)i == e) out.push_back(pe);
    else if ((int)i == sp) out.push_back(ps);
    else out.push_back(pt[i]);
  }
  return true;
}

/*
  RR_MODE 1: the beams of pt rebuilt (on shell, back to back along their own
  axis in the rest frame of the outgoing momenta) so that pt balances
  exactly. The event's momenta balance to ~1e-9 GeV, which is the energy
  scale of its softest photons (x ~ 1e-7); an imbalance of that size moved
  R_2 of such a pair by factors of 2 against S S B.
*/
static bool BalanceOnBeams(Vec4D_Vector &pt, const ATOOLS::Flavour_Vector &fl)
{
  if (pt.size() < 3) return false;
  Vec4D T;
  for (size_t i(2); i < pt.size(); ++i) T += pt[i];
  const double T2(T.Abs2()), m1(fl[0].Mass()), m2(fl[1].Mass());
  if (!(T2 > sqr(m1 + m2)) || !(T[0] > 0.)) return false;
  Poincare rest(T);
  Vec4D a(pt[0]);
  rest.Boost(a);
  Vec3D n(a);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(T2, m1*m1, m2*m2)/T2));
  Vec4D pa(sqrt(m1*m1 + pcm*pcm),  pcm*n), pb(sqrt(m2*m2 + pcm*pcm), -pcm*n);
  rest.BoostBack(pa);
  rest.BoostBack(pb);
  pt[0] = pa;
  pt[1] = pb;
  return true;
}

/*
  RR_MODE 1: the point pt with the initial-state photons in `drop` (indices
  into pt) taken out the way the single real's scaled point does it
  (REAL_MAP 2, MapMomentaScaled): the final system F (final legs and the
  final-state photons, isfsr[slot] for pt[nlegs + slot]) keeps its mass s',
  the kept initial-state photons K and the beams P are scaled by
  x = sqrt(s'/(P - K)^2) and F is boosted rigidly onto x (P - K). Angles to
  the beams and energy fractions of the kept photons, and every invariant
  inside F, are unchanged; with two photons left it is the single real's
  point of the other. Used for the photons other than the pair
  (EventPairPoint): taking a hard one out of its beam instead left that beam
  at a fraction of its energy, and a hard collinear companion inside the
  widened dead cone (4.8e-6 rad, m/E from 1.1e-5 to 5.5e-5) had 400 times
  less crude than in its single real: one pair +171 B in 50k Z-pole nu nu
  events.
*/
static bool ScaleOutISR(const Vec4D_Vector &pt, size_t nlegs, const std::vector<bool> &isfsr,
                        const std::vector<size_t> &drop, const ATOOLS::Flavour_Vector &fl,
                        Vec4D_Vector &out)
{
  if (nlegs < 3 || pt.size() < nlegs || isfsr.size() != pt.size() - nlegs) return false;
  auto dropped = [&](size_t l) {
    return std::find(drop.begin(), drop.end(), l) != drop.end(); };
  for (size_t l : drop) if (l < nlegs || l >= pt.size() || isfsr[l - nlegs]) return false;
  const Vec4D P(pt[0] + pt[1]);
  Vec4D F, K;
  for (size_t l(2); l < pt.size(); ++l) {
    if (dropped(l)) continue;
    if (l < nlegs || isfsr[l - nlegs]) F += pt[l];
    else K += pt[l];
  }
  const Vec4D R(P - K);
  const double sp(F.Abs2()), R2(R.Abs2());
  if (!(sp > 0.) || !(R2 > 0.) || !(R[0] > 0.) || !(F[0] > 0.)) return false;
  const double x(Min(1., sqrt(sp/R2)));
  const double sj(x*x*P.Abs2()), m1(fl[0].Mass()), m2(fl[1].Mass());
  if (!(x > 0.) || !(sj > sqr(m1 + m2))) return false;
  Poincare toP(P);
  Vec4D a(pt[0]);
  toP.Boost(a);
  Vec3D n(a);
  if (!(n.Abs() > 0.)) return false;
  n = n/n.Abs();
  const double pcm(0.5*sqrt(Lambda(sj, m1*m1, m2*m2)/sj));
  Vec4D pa(sqrt(pcm*pcm + m1*m1), pcm*n), pb(sqrt(pcm*pcm + m2*m2), -pcm*n);
  toP.BoostBack(pa);
  toP.BoostBack(pb);
  const Vec4D Fp(pa + pb - x*K);
  if (!(Fp[0] > 0.) || !(Fp.Abs2() > 0.)) return false;
  Poincare fromF(F), toFp(Fp);
  out.assign({pa, pb});
  for (size_t l(2); l < pt.size(); ++l) {
    if (dropped(l)) continue;
    Vec4D q(pt[l]);
    if (l < nlegs || isfsr[l - nlegs]) { fromF.Boost(q); toFp.BoostBack(q); }
    else q = x*q;
    out.push_back(q);
  }
  return BalanceOnBeams(out, fl);
}

/*
  RR_MODE 1: the (n+2)-body point of the pair (i, j) taken from the EVENT:
  the full beams, the post-emission final legs and every photon, with the
  photons other than i and j reduced away: each final-state one recombined
  with its emitter, then the initial-state ones scaled out together
  (ScaleOutISR, the single real's REAL_MAP 2 convention). For a two-photon
  event it is the event itself, so the photons sit at the angles, relative
  to the legs, the generator produced them at (the two-photon map rebuilt
  the pair and put hard collinear photons deep into the dead cone of legs
  they were not generated on). False when the event does not balance.
*/
bool NLO_Base::EventPairPoint(size_t i, size_t j, Vec4D_Vector &P)
{
  const size_t nl(m_flavs.size());
  if (m_postlab.size() != nl || m_bornMomenta.size() < 2) return false;
  Vec4D_Vector ev(m_bornMomenta.begin(), m_bornMomenta.begin() + 2);
  for (size_t l(2); l < nl; ++l) ev.push_back(m_postlab[l]);
  std::vector<size_t> ids;
  for (size_t g(0); g < m_photons.size(); ++g) { ev.push_back(m_photons[g].K()); ids.push_back(g); }
  /*
    The event balances to ~1e-6 GeV (11% of Z-pole mu mu pairs fail
    CheckMomentumConservation's tolerance: p_T up to 4e-6 GeV); the beams are
    rebuilt on their axis in the rest frame of the outgoing momenta, a tilt
    of ~1e-7 rad, far inside every dead cone (m_e/E ~ 1e-5). Larger
    imbalances are not the event.
  */
  {
    Vec4D bal(ev[0] + ev[1]);
    for (size_t l(2); l < ev.size(); ++l) bal -= ev[l];
    if (!(Max(Max(std::abs(bal[0]), std::abs(bal[1])), Max(std::abs(bal[2]), std::abs(bal[3])))
          < 1e-6*sqrt(m_s)) || !BalanceOnBeams(ev, m_flavs)) return false;
  }
  // the other final-state photons recombined with their emitters (last
  // first so that indices stay valid), then the other initial-state photons
  // scaled out together
  for (size_t a(ids.size()); a-- > 0; ) {
    const YFS::Photon &g(m_photons[ids[a]]);
    if (ids[a] == i || ids[a] == j || !g.IsFSR()) continue;
    Vec4D_Vector out;
    bool ok(g.Dip() != nullptr
            && RecombineWithEmitter(ev, nl + a, g.Dip()->Left(), g.Dip()->Right(), m_flavs, out));
    if (!ok) ok = AbsorbPhoton(ev, nl, nl + a, out);
    if (!ok) return false;
    ev = out;
    ids.erase(ids.begin() + a);
  }
  std::vector<bool> isfsr;
  std::vector<size_t> drop;
  for (size_t a(0); a < ids.size(); ++a) {
    isfsr.push_back(m_photons[ids[a]].IsFSR());
    if (ids[a] != i && ids[a] != j) drop.push_back(nl + a);
  }
  if (!drop.empty()) {
    Vec4D_Vector out;
    if (!ScaleOutISR(ev, nl, isfsr, drop, m_flavs, out)) return false;
    ev = out;
    std::vector<size_t> keep;
    for (size_t id : ids) if (id == i || id == j) keep.push_back(id);
    ids = keep;
  }
  if (ev.size() != nl + 2 || ids.size() != 2) return false;
  P = ev;
  if (ids[0] != i) std::swap(P[nl], P[nl + 1]);   // photon i first, j last
  return true;
}

/*
  RR_MODE 1: the pair (i, j)'s exact remainder, see RRMode. P is the event's
  pair point (EventPairPoint; the two-photon MapMomenta point when the event
  does not balance); R_2 is an RR_Generator tree call.
*/
double NLO_Base::RealRealRemainder(size_t i, size_t j)
{
  if (i >= m_rrphot.size() || j >= m_rrphot.size()) return 0.;
  const RRPhotonInfo &gi(m_rrphot[i]), &gj(m_rrphot[j]);
  if (!gi.ok || !gj.ok || !(gi.denom > 0.) || !(gj.denom > 0.) || m_born == 0.) return 0.;
  Vec4D ki(m_photons[i].K()), kj(m_photons[j].K());
  const double phemin(PhotonEminNLO());
  if (phemin > 0. && (ki.E() < phemin || kj.E() < phemin)) return 0.;
  if (Min(ki.E(), kj.E()) < m_rr_soft_cut*sqrt(m_s)) {
    m_softRR++;
    return 0.;
  }
  const size_t nl(m_flavs.size());
  Vec4D_Vector P;
  const bool evpt(EventPairPoint(i, j, P));
  if (!evpt) {
    // fall back to the two-photon map
    P = m_plab;
    MapMomenta(P, ki, kj);
    P.push_back(ki);
    P.push_back(kj);
    CheckMasses(P, 2);
  }
  if (!CheckMomentumConservation(P) || !BalanceOnBeams(P, m_flavs)) { m_zeroRR++; return 0.; }
  const double r2(p_realreal->Calc_R(P)/(2.*pow(2.*M_PI, 6))*BornPhotonSym(2));
  if (p_realreal->FailCut()) { m_failcut = true; return 0.; }
  if (IsBad(r2)) { m_zeroRR++; return 0.; }
  /*
    Each photon's density c is that of its own single-real point: for an
    initial-state photon of a multi-photon event the scaled point (REAL_MAP
    2) carries the photon scaled up, the density ~ 1/x^2 with it. At P it is
    carried by the ratio of the incoherent crude eikonals, photon and legs of
    P over photon and legs of the point (1 when P is the point). Carrying it
    on the generation legs (full beams, m_plab) instead made the double real
    +350% of sigma_fid (a collinear FSR photon against legs it was not
    emitted from); the density on the physical legs is the one the one-photon
    identity uses (SUB8_FSR_LEGS).
  */
  auto carried = [&](const RRPhotonInfo &g, const Vec4D &kP) {
    if (g.pt.size() != m_flavs.size() + 1) return g.denom;
    const double a(CrudeOnLegs(P, kP)), b(CrudeOnLegs(g.pt, g.pt.back()));
    return (a > 0. && b > 0.) ? g.denom*a/b : g.denom; };
  const double ci(carried(gi, P[nl])), cj(carried(gj, P[nl + 1]));
  /*
    The pair's flux: the singles' (Q - k)^2/Q^2 for an initial-state photon
    (full beams Q, FSR photons 1 with REAL_FSR_FLUX 1) taken with both
    initial-state photons of the pair, (Q - K_ISR)^2/Q^2, the event's s'/s
    for two photons. Other flux conventions keep F_i F_j.
  */
  double Fij(gi.flux*gj.flux);
  {
    static const int fsrflux(ATOOLS::Settings::GetMainSettings()["YFS"]
                             ["REAL_FSR_FLUX"].Get<int>());
    if (m_flux_mode == 0 && fsrflux == 1) {
      const Vec4D Q(P[0] + P[1]);
      Vec4D K;
      if (!m_photons[i].IsFSR()) K += P[nl];
      if (!m_photons[j].IsFSR()) K += P[nl + 1];
      if (Q.Abs2() > 0.) Fij = (Q - K).Abs2()/Q.Abs2();
    }
  }
  static const int a0(ATOOLS::Settings::GetMainSettings()["YFS"]
                      ["REAL_ALPHA0"].SetDefault(1).Get<int>());
  const double kap(m_rescale_alpha);
  const double rhoij((a0 ? kap*kap : 1.)*r2*Fij/(ci*cj*m_born));
  double others(1.);
  for (size_t k(0); k < m_rrphot.size(); ++k)
    if (k != i && k != j) others *= m_rrphot[k].factor;
  const double rem(m_born*(rhoij - gi.factor*gj.factor)*others);
  if (m_photons.size() == 2) m_rrExactPair = rhoij;
  static const int probe(ATOOLS::Settings::GetMainSettings()["YFS"]["RR_PROBE"]
                         .SetDefault(0).Get<int>());
  if (probe) {
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ RRPAIR xi="<<2.*m_photons[i].E()/sqrt(m_s)
     <<" xj="<<2.*m_photons[j].E()/sqrt(m_s)<<" fi="<<(m_photons[i].IsFSR()?1:0)
     <<" fj="<<(m_photons[j].IsFSR()?1:0)<<" nph="<<m_photons.size()
     <<" rem/B="<<rem/m_born<<" rhoij="<<rhoij<<" facti="<<gi.factor
     <<" factj="<<gj.factor<<" others="<<others<<" ci/c="<<ci/gi.denom
     <<" cj/c="<<cj/gj.denom<<" Fij/FiFj="<<Fij/(gi.flux*gj.flux)
     <<" evpt="<<(evpt?1:0);
    if (probe >= 2) {
      o<<std::setprecision(17)<<" P=";
      for (const Vec4D &q : P) o<<q<<" ";
      o<<" PB=";
      for (const Vec4D &q : m_plab) o<<q<<" ";
    }
    o<<"\n";
    std::cerr<<o.str();
  }
  if (IsBad(rem)) { m_zeroRR++; return 0.; }
  if (!IsZero(rem)) m_nonZeroRR++;
  return rem;
}

double NLO_Base::CalculateRealReal() {
  m_rrExactPair = std::numeric_limits<double>::quiet_NaN();
  if (!m_rrtool)
    return 0;
  if (RRMode() == 1) {
    static bool checked(false);
    if (!checked) {
      checked = true;
      if (!(m_rr_soft_cut > 0.))
        THROW(fatal_error, "YFS: RR_MODE 1 needs RR_SOFT_CUT > 0 (the lower edge of "
              "the soft x hard form-factor difference in the real-virtual).");
      if (!m_realvirt || RVMode() != 1)
        msg_Error()<<"YFS: RR_MODE 1 without RV_MODE 1: the double real's soft x "
                   <<"hard log is not cancelled; the result depends on RR_SOFT_CUT."
                   <<std::endl;
    }
  }
  double rr(0);
  m_rr_hard2 = 0.;
  const YFS::Photon_Vector &photons(m_photons);
  if (photons.size() == 0)
    return 0;
  const YFS::Photon *h1 = nullptr, *h2 = nullptr;
  for (const YFS::Photon &g : photons) {
    if (!h1 || g.E() > h1->E())      { h2 = h1; h1 = &g; }
    else if (!h2 || g.E() > h2->E()) { h2 = &g; }
  }
  if (RRMode() == 1) {
    for (size_t i(0); i < photons.size(); ++i)
      for (size_t j(i + 1); j < photons.size(); ++j) {
        const double contrib(RealRealRemainder(i, j));
        rr += contrib;
        const YFS::Photon *pi(&photons[i]), *pj(&photons[j]);
        if ((pi == h1 && pj == h2) || (pi == h2 && pj == h1)) m_rr_hard2 = contrib;
      }
    return rr;
  }
  for (int i = 0; i < photons.size(); ++i) {
    for (int j = i + 1; j < photons.size(); ++j) {
      Vec4D k = photons[i].K();
      Vec4D kk = photons[j].K();
      const double phemin(PhotonEminNLO());
      if (phemin>0.0 && (k.E()<phemin || kk.E()<phemin)) continue;
      double contrib = CalculateRealReal(k, kk);
      static const bool betacheck(ATOOLS::Settings::GetMainSettings()["YFS"]["BETA_RECURSION"].Get<int>()!=0);
      if (betacheck) {
        const double gen(CalculateRealN((1u<<i) | (1u<<j)));
        std::cerr<<"@@@ BETA2 hand="<<contrib<<" rec="<<gen
                 <<" dev="<<(contrib!=0. ? std::abs(gen/contrib-1.) : -1.)
                 <<std::endl;
      }
      rr += contrib;
      const YFS::Photon *pi = &photons[i], *pj = &photons[j];
      if ((pi == h1 && pj == h2) || (pi == h2 && pj == h1))
        m_rr_hard2 = contrib;
      if (m_check_rr_sub == 2) {
        // accumulating scatter: record each photon of the pair with the pair
        // residual, so energetic collinear photons in large residuals show up
        RecordSubScatter(k,  contrib, "rr", m_rr_eik);
        RecordSubScatter(kk, contrib, "rr", m_rr_eik);
      }
      if (m_check_rr_sub == 1) {
        // k*=2;
        // kk*=2;
        if (k.E() < 0.2 * sqrt(m_s))
          continue;
        if (kk.E() < 0.2 * sqrt(m_s))
          continue;
        if (!m_failcut)
          CheckRealRealSub(k, kk);
      }
    }
  }
  return rr;
}

double NLO_Base::CalculateRealReal(Vec4D k1, Vec4D k2) {
  if (m_rr_soft_cut > 0.) {
    const double emin = Min(k1.E(), k2.E());
    if (emin < m_rr_soft_cut * sqrt(m_s)) {
      m_softRR++;
      msg_Debugging() << METHOD << " skipping soft photon pair: min(E1,E2)="
                      << emin << " < RR_SOFT_CUT*sqrt(s)="
                      << m_rr_soft_cut * sqrt(m_s) << "\n";
      return 0;
    }
  }
  const double norm = 2. * pow(2 * M_PI, 6);
  m_rr_eik = 0.;  // reset; set once the crude eikonals are computed below
  Vec4D_Vector p(m_plab);
  Vec4D_Vector pp = p;
  Vec4D kk1 = k1, kk2 = k2;

  msg_Debugging() << METHOD
                  << " k1=" << k1 << " E1=" << k1.E() << " pt1=" << k1.PPerp()
                  << " k2=" << k2 << " E2=" << k2.E() << " pt2=" << k2.PPerp() << "\n";

  MapMomenta(p, k1, k2);

  p.push_back(k1);
  p.push_back(k2);

  Vec4D_Vector _p = p;
  _p.pop_back();
  _p.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, _p, m_plab);

  const double subloc1 = p_nlodipoles->CalculateRealSub(k1);
  const double subloc2 = p_nlodipoles->CalculateRealSub(k2);

  double flux;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k1 + k2);
  else
    flux = p_dipoles->CalculateFlux(k1) * p_dipoles->CalculateFlux(k2);

  msg_Debugging() << METHOD << " subloc1=" << subloc1 << " subloc2=" << subloc2
                  << " flux=" << flux << " (mode=" << m_flux_mode << ")\n";

  if (!CheckMomentumConservation(p)) {
    msg_Debugging() << METHOD << " momentum conservation failed, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  double r = p_realreal->Calc_R(p) / norm * BornPhotonSym(2);
  if (p_realreal->FailCut()) {
    msg_Debugging() << METHOD << " FailCut triggered, returning 0\n";
    m_failcut = true;
    m_zeroRR++;
    return 0;
  }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << r << " flux=" << flux << "\n";
    m_zeroRR++;
    return 0;
  }

  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);


  const double sub1  = p_dipoles->CalculateRealSubEEX(kk1);
  const double sub2  = p_dipoles->CalculateRealSubEEX(kk2);
  m_rr_eik = sub1 * sub2;
  const double real1 = CalculateReal(kk1, /*raw*/true);
  const double F1(m_lastflux), SL1(m_lastsubloc), SB1(m_lastdenom);
  const double SL01(m_lastsubloc0), SB01(m_lastdenom0);
  // the eikonal of k2 on the legs of photon 1's one-photon point (the
  // dipoles CalculateReal just built there), and vice versa below
  const double A2raw(p_nlodipoles->CalculateRealSub(kk2));
  const double real2 = CalculateReal(kk2, /*raw*/true);
  const double F2(m_lastflux), SL2(m_lastsubloc), SB2(m_lastdenom);
  const double SL02(m_lastsubloc0), SB02(m_lastdenom0);
  const double A1raw(p_nlodipoles->CalculateRealSub(kk1));
  // what REAL_SUB_EIK / REAL_FSR_FLUX did to each photon's eikonal and
  // denominator in the single real, as a factor (exactly 1 in mode 0)
  auto ratio = [](double a, double b) { return (b != 0. && !IsBad(a/b)) ? a/b : 1.; };
  const double rs1(ratio(SL1, SL01)), rs2(ratio(SL2, SL02));
  const double rd1(ratio(SB1, SB01)), rd2(ratio(SB2, SB02));
  m_recola_evts += 1;
  /*
    YFS: RR_CONVENTIONS (1). beta_2 = R_2 F_1 F_2 - S_2 beta_1(k_1)
    - S_1 beta_1(k_2) - S_1 S_2 B over the crude of each photon, with the flux,
    subtraction eikonal and crude of each photon taken from the single real's
    own calculation of that photon (m_last*). The old assembly (0) built its
    own: the plain eikonal at the double-mapped point, an initial-state flux
    product for every photon and the EEX crude, while beta_1 inside it came
    from the single real. With REAL_SUB_EIK 8 and REAL_FSR_FLUX 1 the two no
    longer cancelled for hard photons: Z-pole mu mu YFS.NLO+RR -187 +- 175%
    and 404 +- 57% pb against 1045 +- 5.7% with the old single-real defaults
    (2026-09-26). The coupling powers are also per term: an eikonal is
    alpha(0), the reals the model's, so S beta_1 carries one m_rescale_alpha
    and S S B two (the old code divided all of it by one).
  */
  /*
    Status 2026-09-26 (Z-pole mu mu, 20k events, Check_RR_Sub 1 for the soft
    limits): 0 is the old assembly. It is small (+0.5%) and passes both
    single-soft limits only with the old single-real settings AND
    REAL_FSR_MAP 0; with REAL_FSR_MAP 2 or 4 (the two-photon FSR map does not
    reduce to the one-photon one) |beta_2|/B stalls at 0.7 as one photon goes
    soft, and the correction is -11%. 3 (default) passes every soft limit with
    REAL_FSR_MAP 0 under REAL_SUB_EIK 8 / REAL_FSR_FLUX 1, but is still -2.8%;
    7 adds the denominator factor (+3.6 +- 3.5%); 8 (eikonals on the other
    photon's legs) made it worse. Open: the FSR maps, and REAL_SUB_EIK 8's
    multichannel crude inside beta_2.
  */
  static const int rrconv(ATOOLS::Settings::GetMainSettings()["YFS"]
                          ["RR_CONVENTIONS"].SetDefault(3).Get<int>());
  if (rrconv) {
    // bitmask while the right combination is established: 1 = the single
    // real's flux, 2 = its subtraction eikonal, 4 = its denominator; a clear
    // bit keeps the old double-real piece (flux product, eikonal at the
    // double-mapped point, EEX crude)
    // 1: the single real's flux; 2: the double real's own eikonal (at its
    // double-mapped point) times the single real's modification factor;
    // 4: the same for the denominator (the EEX crude of the event)
    const double fF(rrconv & 1 ? F1 * F2 : flux);
    const double s1(rrconv & 2 ? subloc1 * rs1 : subloc1), s2(rrconv & 2 ? subloc2 * rs2 : subloc2);
    const double d1(rrconv & 4 ? sub1 * rd1 : sub1), d2(rrconv & 4 ? sub2 * rd2 : sub2);
    if (IsZero(real1) || IsZero(real2) || !(d1 > 0.) || !(d2 > 0.)) {
      m_zeroRR++;
      return 0;
    }
    const double ra(m_rescale_alpha);
    double num(r * fF - (s2 * real1 + s1 * real2) / ra
               - s1 * s2 * m_born / (ra * ra));
    /*
      8: the soft limits. As k2 -> soft, R_2 -> S(k2; legs of photon 1's
      point) R_1(k1), so the eikonal multiplying beta_1(k1) must be k2's on
      THOSE legs (a2), not at the double-mapped point where k1 is removed;
      for FSR the final legs recoil against a hard k1 and the two differ at
      O(1), leaving a finite residual for every soft photon (Z-pole mu mu
      YFS.NLO+RR -11%, nu nu with initial-state radiation only +0.005%).
      beta_2 vanishes in both single-soft limits with
      c = a2 S_1 + a1 S_2 - S_1 S_2, S_i the single real's subtraction.
    */
    if (rrconv & 8) {
      const double a1(A1raw * rs1), a2(A2raw * rs2);
      num = r * fF - (a2 * real1 + a1 * real2) / ra
            - (a2 * SL1 + a1 * SL2 - SL1 * SL2) * m_born / (ra * ra);
    }
    const double tot1(num / (d1 * d2));
    m_rr_eik = d1 * d2;
    if (IsBad(tot1)) {
      msg_Error() << METHOD << " NNLO RR is NaN: r=" << r << " F1=" << F1 << " F2=" << F2
                  << " SL1=" << SL1 << " SL2=" << SL2 << " SB1=" << SB1 << " SB2=" << SB2 << "\n";
      return 0;
    }
    if (!IsZero(tot1)) m_nonZeroRR++;
    return tot1;
  }

  msg_Debugging() << METHOD << " r=" << r
                  << " sub1=" << sub1 << " sub2=" << sub2
                  << " real1=" << real1 << " real2=" << real2
                  << " born=" << m_born << "\n";

  if (IsZero(real1) || IsZero(real2)) {
    msg_Debugging() << METHOD << " real1 or real2 is zero, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  const double fullsub = -subloc2 * real1 - subloc1 * real2 - subloc1 * subloc2 * m_born;
  const double tot     = (r * flux + fullsub / m_rescale_alpha) / sub1 / sub2;

  msg_Debugging() << METHOD << " fullsub=" << fullsub
                  << " r*flux=" << r * flux
                  << " tot=" << tot << "\n";

  if (IsBad(tot))
    msg_Error() << METHOD << " NNLO RR is NaN: r=" << r << " flux=" << flux
                << " fullsub=" << fullsub << " sub1=" << sub1 << " sub2=" << sub2 << "\n";

  if (!IsZero(tot)) m_nonZeroRR++;
  return tot;
}
