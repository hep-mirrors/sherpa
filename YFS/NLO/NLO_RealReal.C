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
  RR_MODE exact / singles_product: the final-state pair (l, r) of ev rebuilt
  as the single real's point does (MapMomentaFSRDipole): at Q_D - K, Q_D the
  pair before it radiated (m_plab), K the kept photons, back to back along the
  event's decay axis. The dipole's other photons are absorbed by this.
*/
static bool RebuildPairAlongDecayAxis(int l, int r, const Vec4D &K,
                                      const Vec4D_Vector &plab, const Vec4D_Vector &postlab,
                                      const ATOOLS::Flavour_Vector &fl, Vec4D_Vector &ev)
{
  if (l < 2 || r < 2 || l == r || l >= (int)fl.size() || r >= (int)fl.size()
      || plab.size() != fl.size() || postlab.size() != fl.size() || ev.size() < fl.size())
    return false;
  const Vec4D Qd(plab[l] + plab[r]);
  const Vec4D Qp(Qd - K);
  const double Mp2(Qp.Abs2());
  if (!(Mp2 > 0.) || !(Qp[0] > 0.)) return false;
  const double Mp(sqrt(Mp2));
  const double m1(fl[l].Mass()), m2(fl[r].Mass());
  if (m1 + m2 >= Mp) return false;
  Vec4D q1(postlab[l]), q2(postlab[r]);
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
  ev[l] = f1; ev[r] = f2;
  return true;
}

/*
  RR_MODE 1: the beams of pt rebuilt on shell, back to back along their own
  axis in the rest frame of the outgoing momenta, so that pt balances exactly
  (the event balances only to ~1e-9 GeV, the scale of its softest photons).
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
  into pt) taken out as the single real's scaled point does (MapMomentaScaled):
  the final system F (isfsr[slot] for pt[nlegs + slot]) keeps its mass s', the
  kept ISR photons K and beams P are scaled by x = sqrt(s'/(P - K)^2) and F
  is boosted onto x (P - K). Angles, energy fractions and invariants inside F
  are unchanged. Used for the photons other than the pair (EventPairPoint).
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
  RR_MODE exact / singles_product: the (n+2)-body point of the pair (i, j) on
  the single reals' chain. Other final-state photons are taken out by
  rebuilding their dipole's pair at Q_D minus the pair's own photons of that
  dipole along the event's decay axis (RebuildPairAlongDecayAxis), initial-
  state ones by scaling out (ScaleOutISR). As one photon goes soft the point
  tends to the other's single-real point, so rho_ij -> rho_i rho_j. For a
  two-photon event it is the event itself. Recombining the companion with its
  emitter instead would rotate the pair's axis and bias B(P)/B. False when a
  pair photon has no dipole or the reduction fails.
*/
bool NLO_Base::EventPairPoint(size_t i, size_t j, Vec4D_Vector &P)
{
  return EventSetPoint(std::vector<size_t>{i, j}, P);
}

bool NLO_Base::EventSetPoint(const std::vector<size_t> &S, Vec4D_Vector &P)
{
  const size_t nl(m_flavs.size());
  if (m_postlab.size() != nl || m_plab.size() != nl || m_bornMomenta.size() < 2) return false;
  if (p_dipoles && p_dipoles->PoleActive()) return false;
  auto inset = [&](size_t g) { return std::find(S.begin(), S.end(), g) != S.end(); };
  // the generation legs with the full beams
  Vec4D_Vector ev(m_bornMomenta.begin(), m_bornMomenta.begin() + 2);
  for (size_t l(2); l < nl; ++l) ev.push_back(m_plab[l]);
  // the set's final-state photons: their dipole's pair rebuilt without
  // them, which absorbs that dipole's other photons; a dipole with no set
  // photon keeps its generation legs, which absorbs all of its photons
  std::map<std::pair<int,int>, Vec4D> ksel;
  for (size_t g : S) {
    const YFS::Photon &ph(m_photons[g]);
    if (!ph.IsFSR()) continue;
    if (ph.Dip() == nullptr) return false;
    const int l(ph.Dip()->Left()), r(ph.Dip()->Right());
    ksel[std::make_pair(std::min(l, r), std::max(l, r))] += ph.K();
  }
  for (const auto &sel : ksel)
    if (!RebuildPairAlongDecayAxis(sel.first.first, sel.first.second, sel.second,
                                   m_plab, m_postlab, m_flavs, ev)) return false;
  // the set's photons and the initial-state companions, the latter to be
  // scaled out; ordg: the set's photons in m_photons order
  std::vector<bool> isfsr;
  std::vector<size_t> drop, ordg;
  for (size_t g(0); g < m_photons.size(); ++g) {
    const bool mem(inset(g));
    if (!mem && m_photons[g].IsFSR()) continue;
    if (mem) ordg.push_back(g);
    if (!mem) drop.push_back(ev.size());
    isfsr.push_back(m_photons[g].IsFSR());
    ev.push_back(m_photons[g].K());
  }
  Vec4D_Vector out;
  if (!ScaleOutISR(ev, nl, isfsr, drop, m_flavs, out)) return false;
  if (out.size() != nl + S.size() || ordg.size() != S.size()) return false;
  // out carries the set's photons in m_photons order; P in S order
  P.assign(out.begin(), out.begin() + nl);
  for (size_t g : S)
    P.push_back(out[nl + (std::find(ordg.begin(), ordg.end(), g) - ordg.begin())]);
  return true;
}

/*
  RR_MODE exact / singles_product: the pair (i, j)'s remainder, see RRMode
  (NLO_Base_Internal.H). P is the event's pair point (EventPairPoint; the
  two-photon MapMomenta point when the event does not balance); R_2 is an
  RR_Generator tree call.
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
    Each photon's density c is that of its own single-real point (for an ISR
    photon of a multi-photon event the scaled point carries the 1/x^2). At P
    it is carried by the ratio of the incoherent crude eikonals, photon and
    legs of P over those of the point. Using the generation legs instead would
    put a collinear FSR photon against legs it was not emitted from
    (SUB8_FSR_LEGS).
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
  if (m_flux_mode == fluxmode::event && RealFSRFlux() == realfsrflux::no_flux) {
    const Vec4D Q(P[0] + P[1]);
    Vec4D K;
    if (!m_photons[i].IsFSR()) K += P[nl];
    if (!m_photons[j].IsFSR()) K += P[nl + 1];
    if (Q.Abs2() > 0.) Fij = (Q - K).Abs2()/Q.Abs2();
  }
  const double kap(m_rescale_alpha);
  const double rhoij((RealAlpha0() ? kap*kap : 1.)*r2*Fij/(ci*cj*m_born));
  double others(1.);
  for (size_t k(0); k < m_rrphot.size(); ++k)
    if (k != i && k != j) others *= m_rrphot[k].factor;
  // what the pair's exact over crude is taken beyond: the product's factors
  // (exact) or the singles' own exact over crude (singles_product)
  const bool singles(RRMode() == rrmode::singles_product);
  const double rem(singles ? m_born*(rhoij - gi.rho*gj.rho)*others
                           : m_born*(rhoij - gi.factor*gj.factor)*others);
  if (IsBad(rem)) { m_zeroRR++; return 0.; }
  if (!IsZero(rem)) m_nonZeroRR++;
  return rem;
}

bool NLO_Base::RRSinglesProduct() const
{
  return m_rrtool && RRMode() == rrmode::singles_product;
}

double NLO_Base::CalculateRealReal() {
  if (!m_rrtool)
    return 0;
  if (RRMode() == rrmode::exact) {
    static bool checked(false);
    if (!checked) {
      checked = true;
      if (!(m_rr_soft_cut > 0.))
        THROW(fatal_error, "YFS: RR_MODE exact needs RR_SOFT_CUT > 0 (the lower edge of "
              "the soft x hard form-factor difference in the real-virtual).");
      if (!m_realvirt || !RVIsRemainder())
        msg_Error()<<"YFS: RR_MODE exact without RV_MODE remainder: the double real's soft x "
                   <<"hard log is not cancelled; the result depends on RR_SOFT_CUT."
                   <<std::endl;
    }
  }
  if (RRMode() == rrmode::singles_product) {
    static bool checked(false);
    if (!checked) {
      checked = true;
      if (m_rr_soft_cut > 0.)
        msg_Info()<<"YFS: RR_MODE singles_product with RR_SOFT_CUT "<<m_rr_soft_cut
                  <<": pairs below it keep the product's factorised double real; "
                  <<"0 takes every generated pair."<<std::endl;
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
  if (RRFromSingles()) {
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
  if (m_flux_mode == fluxmode::mapped)
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
    YFS: RR_CONVENTIONS (default single_flux, single_subtraction), a list of
    names (or an integer mask, the value in brackets), each taking one piece
    from the single real's own calculation of each photon (m_last*). beta_2 =
    R_2 F_1 F_2 - S_2 beta_1(k_1) - S_1 beta_1(k_2) - S_1 S_2 B over the crude
    of each photon. A piece not named is the double real's own (flux product,
    eikonal at the double-mapped point, EEX crude).
      single_flux (1): the single real's flux;
      single_subtraction (2): the double real's eikonal times the single real's
        modification factor;
      single_denominator (4): the same for the denominator;
      soft_limit_legs (8): the soft-limit eikonals, below.
    The coupling powers are per term: an eikonal is alpha(0), the reals the
    model's, so S beta_1 carries one m_rescale_alpha and S S B two.
  */
  static const rrconventions rrconv(ReadRRConventions(
      ATOOLS::Settings::GetMainSettings()["YFS"]["RR_CONVENTIONS"]
      .SetDefault(std::vector<std::string>{"single_flux", "single_subtraction"})
      .GetVector<std::string>()));
  if (!rrconv.Legacy()) {
    const double fF(rrconv.single_flux ? F1 * F2 : flux);
    const double s1(rrconv.single_subtraction ? subloc1 * rs1 : subloc1);
    const double s2(rrconv.single_subtraction ? subloc2 * rs2 : subloc2);
    const double d1(rrconv.single_denominator ? sub1 * rd1 : sub1);
    const double d2(rrconv.single_denominator ? sub2 * rd2 : sub2);
    if (IsZero(real1) || IsZero(real2) || !(d1 > 0.) || !(d2 > 0.)) {
      m_zeroRR++;
      return 0;
    }
    const double ra(m_rescale_alpha);
    double num(r * fF - (s2 * real1 + s1 * real2) / ra
               - s1 * s2 * m_born / (ra * ra));
    /*
      soft_limit_legs (8): as k2 -> soft, R_2 -> S(k2; legs of photon 1's
      point) R_1(k1), so the eikonal multiplying beta_1(k1) must be k2's on
      THOSE legs (a2); for FSR the legs recoil against a hard k1 and the two
      differ at O(1). beta_2 then vanishes in both single-soft limits with
      c = a2 S_1 + a1 S_2 - S_1 S_2, S_i the single real's subtraction.
    */
    if (rrconv.soft_limit_legs) {
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
