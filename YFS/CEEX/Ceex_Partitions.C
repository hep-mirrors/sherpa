#include "YFS/CEEX/Ceex_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

#include <algorithm>

using namespace YFS;


bool Ceex_Base::BuildStages() {
  /*
    Stages are the charge-neutral groups of legs that radiate coherently and
    whose propagators the partition sum moves: the initial state, and one
    group per resonance in the final state. The final groups come from YFS's
    radiating dipoles (SetStageGroups); charged final legs in no group form
    one extra group; the initial state is appended LAST, so with a single
    final group this is leg for leg the old (final = 0, initial = 1) split
    and the partition enumeration order is unchanged. With no groups given
    every charged final leg goes into one flat final stage, which is valid
    for any N but resums initial-final interference through the production
    line only.

    theta = -1 incoming, +1 outgoing, so w = Q*theta reproduces the weights
    Sfactor applies implicitly.
  */
  if (m_flavs.size() < 4) return false;
  m_stagemom.clear();
  m_resonances.clear();
  for (size_t i(0); i < m_flavs.size(); ++i)
    m_stagemom.push_back(i < 2 ? m_bornmomenta[i] : m_pceex[i]);

  /*
    The handler's groups in one form: the plain leg lists (SetStageGroups
    with leg indices, every process so far) become StageGroups without a
    resonance, so the code below is the same for both.
  */
  std::vector<StageGroup> handed(m_wgroups);
  for (size_t g(0); g < m_groups.size(); ++g) {
    StageGroup sg;
    sg.externalLegs = m_groups[g];
    handed.push_back(sg);
  }

  // one final group per handed group: its charged, unused external legs and
  // - for a decay stage - its resonance, registered as a new m_stagemom entry
  struct FinalGroup { std::vector<int> legs; int resonance = -1; };
  std::vector<FinalGroup> groups;
  std::vector<char> used(m_flavs.size(), 0);
  for (size_t g(0); g < handed.size(); ++g) {
    FinalGroup grp;
    for (size_t k(0); k < handed[g].externalLegs.size(); ++k) {
      const int l(handed[g].externalLegs[k]);
      if (l < 2 || l >= (int)m_flavs.size() || used[l]) continue;
      if (m_flavs[l].Charge() == 0.) continue;
      grp.legs.push_back(l); used[l] = 1;
    }
    if (handed[g].HasResonance()) {
      Resonance r;
      r.flav = handed[g].resonance;
      r.daughters = handed[g].resonanceDaughters;
      r.stagemomIndex = (int)m_stagemom.size();
      Vec4D pw;
      for (size_t k(0); k < r.daughters.size(); ++k) pw += m_pceex[r.daughters[k]];
      m_stagemom.push_back(pw);
      m_resonances.push_back(r);
      grp.resonance = (int)m_resonances.size() - 1;
    }
    if (!grp.legs.empty() || grp.resonance >= 0) groups.push_back(grp);
  }
  FinalGroup rest;
  for (size_t l(2); l < m_flavs.size(); ++l)
    if (!used[l] && m_flavs[l].Charge() != 0.) rest.legs.push_back((int)l);
  if (!rest.legs.empty()) groups.push_back(rest);
  if (groups.empty()) groups.push_back(FinalGroup());   // neutral final state
  /*
    A stage factor is gauge invariant only if the stage is separately
    CHARGE-NEUTRAL. A dipole group that is not (a same-charge pair from the
    selector's fallback pass) is folded back into one flat final stage; if
    even that is charged the decomposition is refused, as before. A decay
    stage counts its resonance as an INCOMING leg: {mu+, W+ in} is neutral.
  */
  auto charge = [&](const FinalGroup &g) {
    double q(0.);
    for (size_t k(0); k < g.legs.size(); ++k) q += m_flavs[g.legs[k]].Charge();
    if (g.resonance >= 0) q -= m_resonances[g.resonance].flav.Charge();
    return q;
  };
  bool neutral(true);
  for (size_t g(0); g < groups.size(); ++g) if (!IsZero(charge(groups[g]))) neutral = false;
  if (!neutral && groups.size() > 1) {
    static bool warned(false);
    if (!warned) { warned = true;
      msg_Debugging()<<METHOD<<"(): a final-state group carries net charge; using"
                <<" one flat final stage instead.\n"; }
    FinalGroup all;
    for (size_t l(2); l < m_flavs.size(); ++l)
      if (m_flavs[l].Charge() != 0.) all.legs.push_back((int)l);
    groups.assign(1, all);
    m_resonances.clear();
    m_stagemom.resize(m_flavs.size());
  }
  m_stagelegs.clear();
  for (size_t g(0); g < groups.size(); ++g) {
    std::vector<StageLeg> st;
    for (size_t k(0); k < groups[g].legs.size(); ++k)
      st.push_back({groups[g].legs[k], m_flavs[groups[g].legs[k]].Charge()});
    if (groups[g].resonance >= 0) {
      Resonance &r(m_resonances[groups[g].resonance]);
      r.decayStage = (int)m_stagelegs.size();
      st.push_back({r.stagemomIndex, -r.flav.Charge()});   // incoming
    }
    m_stagelegs.push_back(st);
  }
  /*
    The production stage: the beams, and every resonance as an OUTGOING
    leg. With no resonance this is the initial stage as before. The W terms
    cancel between here and the decay stages in the total eikonal, so the
    closure below holds with any common W momentum.
  */
  { std::vector<StageLeg> ini;
    for (int i(0); i < 2; ++i) ini.push_back({i, -m_flavs[i].Charge()});
    for (size_t r(0); r < m_resonances.size(); ++r)
      ini.push_back({m_resonances[r].stagemomIndex, m_resonances[r].flav.Charge()});
    m_stagelegs.push_back(ini); }
  m_initstage = (int)m_stagelegs.size() - 1;
  for (size_t g(0); g < m_stagelegs.size(); ++g) {
    double qsum(0.);
    for (size_t i(0); i < m_stagelegs[g].size(); ++i) qsum += m_stagelegs[g][i].w;
    if (!IsZero(qsum)) {
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): stage "<<g<<" carries net charge "<<qsum
                   <<", so its soft factor is gauge dependent. A valid stage "
                   <<"decomposition partitions the charged legs into "
                   <<"charge-neutral subsets. Refusing this decomposition."
                   <<std::endl;
      }
      return false;
    }
  }
  m_nstages = (int)m_stagelegs.size();
  m_stagereduces.assign(m_nstages, 0);
  m_stagereduces[m_initstage] = 1;   // only initial-stage photons reduce X
  { static bool shown(false);
    if (!shown) { shown = true;
      msg_Debugging()<<METHOD<<"(): "<<m_nstages<<" stages:";
      for (size_t g(0); g < m_stagelegs.size(); ++g) {
        msg_Debugging()<<" {";
        for (size_t i(0); i < m_stagelegs[g].size(); ++i) {
          const int leg(m_stagelegs[g][i].leg);
          if (i) msg_Debugging()<<",";
          if (IsResonanceLeg(leg))
            msg_Debugging()<<ResonanceOf(leg).flav
                           <<(m_stagelegs[g][i].w*ResonanceOf(leg).flav.Charge() < 0.
                              ? "(in)" : "(out)");
          else msg_Debugging()<<leg;
        }
        msg_Debugging()<<"}"<<(g==(size_t)m_initstage?"(production)":""); }
      msg_Debugging()<<"\n"; } }
  if (WStagesActive()) {
    static bool shown(false);
    if (!shown) { shown = true;
      msg_Info()<<"CEEX: W stages active for "<<m_flavs.size()-2
                <<"-fermion final state: "<<m_nstages<<" stages (production with "
                <<m_resonances.size()<<" W legs, one decay stage per W)."
                <<std::endl; }
  }
  return true;
}


bool Ceex_Base::IsResonanceLeg(int leg) const
{
  return leg >= (int)m_flavs.size() && leg < (int)m_stagemom.size();
}

const Ceex_Base::Resonance &Ceex_Base::ResonanceOf(int leg) const
{
  for (size_t r(0); r < m_resonances.size(); ++r)
    if (m_resonances[r].stagemomIndex == leg) return m_resonances[r];
  THROW(fatal_error, "stage leg "+ToString(leg)+" names no resonance");
}

std::vector<int> Ceex_Base::StageExternalLegs(int stage) const
{
  std::vector<int> legs;
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return legs;
  const std::vector<StageLeg> &L(m_stagelegs[stage]);
  // insertion order, each leg once: the momentum sums built from this list
  // then add in the same order as before, which the bit-identity of the
  // no-resonance case depends on
  auto add = [&](int leg) {
    if (std::find(legs.begin(), legs.end(), leg) == legs.end()) legs.push_back(leg); };
  for (size_t l(0); l < L.size(); ++l) {
    if (!IsResonanceLeg(L[l].leg)) { add(L[l].leg); continue; }
    const std::vector<int> &d(ResonanceOf(L[l].leg).daughters);
    for (size_t k(0); k < d.size(); ++k) add(d[k]);
  }
  return legs;
}

size_t Ceex_Base::StageShiftMask(int stage) const
{
  size_t mask(0);
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return mask;
  if (stage == m_initstage) {
    for (size_t l(0); l < m_stagelegs[stage].size(); ++l) {
      const int leg(m_stagelegs[stage][l].leg);
      if (leg >= 0 && leg < 2) mask |= ((size_t)1) << leg;
    }
    return mask;
  }
  const std::vector<int> legs(StageExternalLegs(stage));
  for (size_t l(0); l < legs.size(); ++l)
    if (legs[l] >= 0 && legs[l] < (int)m_flavs.size()) mask |= ((size_t)1) << legs[l];
  return mask;
}

Vec4D Ceex_Base::StageSystemMomentum(int stage) const
{
  Vec4D q;
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return q;
  const std::vector<int> legs(StageExternalLegs(stage));
  for (size_t l(0); l < legs.size(); ++l)
    if (legs[l] >= 0 && legs[l] < (int)m_pceex.size()) q += m_pceex[legs[l]];
  return q;
}

Vec4D Ceex_Base::StagePhotonSum(int stage, int iphot, int lphot) const
{
  Vec4D K;
  if (m_stage.size() != m_allphotons.size()) return K;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot && (int)i != lphot && m_stage[i] == stage)
      K += m_allphotons[i];
  return K;
}

Ceex_Base::LegMomentum
Ceex_Base::StageLegMomentum(const StageLeg &l, int iphot,
                            const Vec4D_Vector &p) const
{
  LegMomentum m;
  if (!IsResonanceLeg(l.leg)) {
    if (l.leg >= 0 && l.leg < (int)p.size()) m.num = m.pole = p[l.leg];
    return m;
  }
  const Resonance &r(ResonanceOf(l.leg));
  for (size_t k(0); k < r.daughters.size(); ++k)
    if (r.daughters[k] < (int)p.size()) m.num += p[r.daughters[k]];
  m.pole = m.num;
  if (m_weikonal != 0 && r.decayStage >= 0)
    m.pole += StagePhotonSum(r.decayStage, iphot);
  return m;
}

Complex Ceex_Base::StageCurrent(int stage, int iphot, const Vec4D &k, int hel,
                                const Vec4D_Vector &p)
{
  Complex tot(0., 0.);
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return tot;
  const std::vector<StageLeg> &L(m_stagelegs[stage]);
  for (size_t l(0); l < L.size(); ++l) {
    if (L[l].w == 0.) continue;
    const LegMomentum m(StageLegMomentum(L[l], iphot, p));
    tot += L[l].w * SfactorLeg(m.num, m.pole, k, hel);
  }
  return tot;
}

void Ceex_Base::RecomputeResonanceSfactors()
{
  if (!WStagesActive() || m_weikonal == 0) return;
  if ((int)m_Sfac.size() != m_nstages) return;
  for (int g(0); g < m_nstages; ++g) {
    bool hasres(false);
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (IsResonanceLeg(m_stagelegs[g][l].leg)) hasres = true;
    if (!hasres) continue;
    for (size_t i(0); i < m_allphotons.size() && i < m_Sfac[g].size(); ++i)
      m_Sfac[g][i] = StageCurrent(g, (int)i, m_allphotons[i], m_PhoHel[i], m_pceex);
  }
}


void Ceex_Base::CalculateSfactors() {
  if (!BuildStages()) { m_Sfac.clear(); return; }
  m_Sfac.assign(m_nstages, std::vector<Complex>());
  for (size_t i(0); i < m_allphotons.size(); ++i)
    for (int g(0); g < m_nstages; ++g) {
      /*
        The stage's eikonal current, summed over its own legs against their
        weights w = Q*theta.
        
        Legs name entries of m_stagemom, which holds the external legs first
        and any reconstructed resonance after them.
      */
      const std::vector<StageLeg> &L(m_stagelegs[g]);
      Complex sg(0., 0.);
      for (size_t l(0); l < L.size(); ++l)
        sg += L[l].w * SfactorLeg(m_stagemom[L[l].leg], m_allphotons[i],
                                  m_PhoHel[i]);
      m_Sfac[g].push_back(sg);
    }
  // a resonance leg sits at its daughters here (no partition yet); under
  // W_EIKONAL 1 the partition loop moves its pole per partition

  /*
    Closure the stage currents must sum to the TOTAL eikonal,

        sum_g s_g(k)  =  s(k)  =  sum_{all legs i} w_i * SfactorLeg(p_i,k)

    which is the statement that every emitter is counted exactly once. 
  */
  if (m_checkxs) {
    for (size_t i(0); i < m_allphotons.size(); ++i) {
      Complex tot(0., 0.), sum(0., 0.);
      for (size_t l(0); l < m_flavs.size() && l < m_stagemom.size(); ++l) {
        const double w(m_flavs[l].Charge() * (l < 2 ? -1. : +1.));
        tot += w * SfactorLeg(m_stagemom[l], m_allphotons[i], m_PhoHel[i]);
      }
      for (int g(0); g < m_nstages; ++g) sum += m_Sfac[g][i];
      const double den(std::abs(tot) + std::abs(sum));
      const double rel(den > 0. ? std::abs(sum - tot)/den : 0.);
      if (rel > 1e-10) {
        static bool warned(false);
        if (!warned) {
          warned = true;
          msg_Error()<<METHOD<<"(): stage closure violated, |sum_g s_g - s| / "
                     <<"scale = "<<rel<<". The stage decomposition is counting "
                     <<"an emitter twice or missing one."<<std::endl;
        }
      }
    }
  }
}


/*
  Exact momentum conservation of what Comix is handed.

  Comix's Berends-Giele recursion builds every internal line from the legs on
  the non-root side. For the beam-electron (root) emission line in M_1 that
  is the COMPLEMENT of {e-, gamma}: the sum of the positron, every final leg
  and - through the partition shifts - every other photon. That sum equals
  p_e - k only if the set conserves momentum exactly, and the collinear
  virtuality it has to reproduce is E_k E (theta^2 + m^2/E^2): 1e-8 GeV^2
  for a 0.2 GeV photon at 25 microradians. The event as handed to CEEX
  balances to 1e-9 GeV at the median and 3e-7 at the 90th percentile
  (four-fermion final states at 250 GeV; the FSR mapping's tolerance, not
  rounding), which is a 1e-6..1e-4 GeV^2 error on that line - larger than
  the line itself. Comix's own ProjectWideMomenta would fix it, but the CEEX
  entry points switch it off, rightly: their leg set is missing the other
  photons on purpose and the projection would move the legs by THEIR
  momenta.

  Measured (four-fermion, 250 GeV, 470 beam-collinear enumerated photons):
  |M_1|/|s B| on the POSITRON side, whose emission line is the direct sum
  {e+, gamma}, is 1.00 at every angle down to 3e-5 rad; on the electron side
  it is 0.99 when the event's residual is below 1e-12 GeV and 0.02 when it is
  1e-8..1e-6, median 0.25 below 3e-5 rad. Each such photon then contributes
  beta_1 = M_1 - s B = -s B and cancels the whole Born: 3.6% of events had
  rho_1/rho_0 < 0.1, 16% a CEEX factor below 0.2. Same in mu mu nu nu and
  in 2 -> 2.

  So the set is projected here, once per event, before the eikonals, the
  partition Borns and M_1 read it: minimise sum |delta_i|^2 (Euclidean, over
  the final legs and every photon; the beams stay as they are, exact and on
  axis) subject to sum delta_i = R and 2 p_i.delta_i = m_i^2 - p_i^2, the
  same problem ProjectWideMomenta solves, as a 4 x 4 linear system for the
  Lagrange vector, two Newton passes (the mass-shell constraint is
  quadratic; after the first pass it is violated by |delta|^2 ~ 1e-14 GeV^2,
  after the second by 1e-28). Corrections are the residual's size, 1e-7 GeV
  at most: nothing physical moves, and a point that already balances to
  rounding is returned unchanged to rounding.
*/
bool Ceex_Base::RepairMomentumBalance()
{
  static const int on(ATOOLS::Settings::GetMainSettings()["CEEX"]
                      ["MOMENTUM_REPAIR"].Get<int>());
  m_repairresid = -1.;
  if (!on || m_pceex.size() < 3) return false;
  std::vector<Vec4D*> v;
  std::vector<double> m2;
  for (size_t i(2); i < m_pceex.size(); ++i) {
    v.push_back(&m_pceex[i]);
    m2.push_back(sqr(i < m_flavs.size() ? m_flavs[i].Mass() : 0.));
  }
  for (size_t i(0); i < m_isrphotons.size(); ++i) { v.push_back(&m_isrphotons[i]); m2.push_back(0.); }
  for (size_t i(0); i < m_fsrphotons.size(); ++i) { v.push_back(&m_fsrphotons[i]); m2.push_back(0.); }
  const size_t n(v.size());
  if (n < 2) return false;
  const Vec4D P(m_pceex[0] + m_pceex[1]);
  const double scale(P[0] > 0. ? P[0] : 1.);
  for (int pass(0); pass < 3; ++pass) {
    Vec4D R(P);
    for (size_t i(0); i < n; ++i) R -= *v[i];
    double rmax(0.);
    std::vector<double> r(n);
    for (size_t i(0); i < n; ++i) {
      r[i] = m2[i] - v[i]->Abs2();
      rmax = Max(rmax, dabs(r[i]));
    }
    const double rn(sqrt(sqr(R[0]) + sqr(R[1]) + sqr(R[2]) + sqr(R[3])));
    if (pass == 0) m_repairresid = rn;
    // already at rounding: leave the point alone
    if (rn < 1e-15*scale && rmax < 1e-15*scale*scale) break;
    // A Lambda = b with A = n 1 - sum_i u_i u_i^T / |u_i|^2, u_i = G p_i
    double A[4][4], b[4];
    for (int a(0); a < 4; ++a) {
      b[a] = R[a];
      for (int c(0); c < 4; ++c) A[a][c] = (a == c ? double(n) : 0.);
    }
    std::vector<double> inv(n);
    for (size_t i(0); i < n; ++i) {
      const Vec4D &p(*v[i]);
      const double u[4] = {p[0], -p[1], -p[2], -p[3]};
      const double nn(u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
      inv[i] = nn > 0. ? 1./nn : 0.;
      for (int a(0); a < 4; ++a) {
        b[a] -= 0.5*r[i]*u[a]*inv[i];
        for (int c(0); c < 4; ++c) A[a][c] -= u[a]*u[c]*inv[i];
      }
    }
    // Gaussian elimination with partial pivoting
    double L[4] = {0., 0., 0., 0.};
    { int piv[4] = {0, 1, 2, 3};
      for (int c(0); c < 4; ++c) {
        int best(c);
        for (int a(c + 1); a < 4; ++a)
          if (dabs(A[piv[a]][c]) > dabs(A[piv[best]][c])) best = a;
        std::swap(piv[c], piv[best]);
        const double d(A[piv[c]][c]);
        if (d == 0.) return false;
        for (int a(c + 1); a < 4; ++a) {
          const double f(A[piv[a]][c]/d);
          if (f == 0.) continue;
          for (int cc(c); cc < 4; ++cc) A[piv[a]][cc] -= f*A[piv[c]][cc];
          b[piv[a]] -= f*b[piv[c]];
        }
      }
      for (int c(3); c >= 0; --c) {
        double x(b[piv[c]]);
        for (int cc(c + 1); cc < 4; ++cc) x -= A[piv[c]][cc]*L[cc];
        L[c] = x/A[piv[c]][c];
      } }
    for (size_t i(0); i < n; ++i) {
      Vec4D &p(*v[i]);
      const double u[4] = {p[0], -p[1], -p[2], -p[3]};
      const double mu((0.5*r[i] - (u[0]*L[0] + u[1]*L[1] + u[2]*L[2] + u[3]*L[3]))*inv[i]);
      for (int a(0); a < 4; ++a) p[a] += L[a] + mu*u[a];
    }
  }
  return true;
}


void Ceex_Base::PartitionStart(int &last) {
  const bool fsr(HasFSR());
  static const double softcut(ATOOLS::Settings::GetMainSettings()["CEEX"]
                              ["SOFT_PARTITION_CUT"].Get<double>());
  const size_t n(m_allphotons.size());
  m_stage.assign(n, fsr ? 0 : m_initstage);
  m_fixedstage.assign(n, 0);
  size_t nfree(0);
  for (size_t i(0); i < n; ++i) {
    const double x(m_s > 0. ? 2.*m_allphotons[i][0]/sqrt(m_s) : 1.);
    if (fsr && x >= softcut && (int)m_Sfac.size() == m_nstages) { ++nfree; continue; }
    // fixed: the stage with the largest eikonal is its natural emitter
    m_fixedstage[i] = 1;
    int best(m_initstage); double bmax(-1.);
    if ((int)m_Sfac.size() == m_nstages)
      for (int g(0); g < m_nstages; ++g)
        if (i < m_Sfac[g].size() && std::abs(m_Sfac[g][i]) > bmax)
          { bmax = std::abs(m_Sfac[g][i]); best = g; }
    m_stage[i] = fsr ? best : m_initstage;
  }
  last = (fsr && nfree > 0) ? 0 : 1;
}


void Ceex_Base::PartitionPlus(int &last) {
  const size_t n(m_stage.size());
  if (n == 0) { last = 2; return; }
  // Odometer in base m_nstages over the ENUMERATED photons only, least
  // significant photon first. At m_nstages = 2 with every photon enumerated
  // this visits partitions in the same order as the original ISR/FSR
  // increment it replaces.
  size_t i(0);
  for (; i < n; ++i) {
    if (m_fixedstage.size() == n && m_fixedstage[i]) continue;
    if (++m_stage[i] < m_nstages) break;
    m_stage[i] = 0;
  }
  if (i == n) { last = 2; return; }   // carried off the top: enumeration done
  // last = 1 marks the final partition, so the caller processes it and stops.
  bool atmax(true);
  for (size_t j(0); j < n; ++j) {
    if (m_fixedstage.size() == n && m_fixedstage[j]) continue;
    if (m_stage[j] != m_nstages - 1) { atmax = false; break; }
  }
  last = atmax ? 1 : 0;
}


void Ceex_Base::PrintKinematics(std::ostream &o) const
{
  o<<std::setprecision(6)<<"HEAVY legs:";
  for (size_t i(0); i < m_pceex.size(); ++i)
    o<<" ["<<i<<"] "<<m_flavs[i]<<" E="<<m_pceex[i][0]
     <<" th="<<m_pceex[i].Theta()<<" m="<<m_pceex[i].Mass();
  o<<std::endl;
  const double rs(m_s > 0. ? sqrt(m_s) : 1.);
  for (size_t j(0); j < m_allphotons.size(); ++j) {
    const Vec4D &k(m_allphotons[j]);
    double thmin(M_PI);
    for (size_t i(2); i < m_pceex.size(); ++i)
      if (m_flavs[i].Charge() != 0.) {
        const double c(Vec3D(k)*Vec3D(m_pceex[i])
                       /(Vec3D(k).Abs()*Vec3D(m_pceex[i]).Abs()));
        thmin = Min(thmin, acos(Max(-1., Min(1., c))));
      }
    o<<"HEAVY photon "<<j<<" E="<<k[0]<<" x="<<2.*k[0]/rs
     <<" th="<<k.Theta()<<" thF="<<thmin
     <<" hel="<<(j < m_PhoHel.size() ? m_PhoHel[j] : 0)
     <<" fixed="<<(j < m_fixedstage.size() ? (int)m_fixedstage[j] : -1);
    for (int g(0); g < m_nstages; ++g)
      if ((int)m_Sfac.size() == m_nstages && j < m_Sfac[g].size())
        o<<" |s"<<g<<"|="<<std::abs(m_Sfac[g][j]);
    o<<std::endl;
  }
  o<<"HEAVY svarQ="<<m_svarQ<<" svarY="<<m_svarY<<" s="<<m_s
   <<" rhocrud="<<m_rhocrud<<" result0="<<m_result0<<std::endl;
  // the two physical t-channel invariants (Bhabha)
  o<<"@@@ HEAVY t02="<<T02()<<" t13="<<T13()<<std::endl;
}


void Ceex_Base::MakePhotonHel() {
  // CEEX_PIN_PHOTON_HEL pins the helicity for validation. The random number is
  // still drawn, so the RNG stream - and therefore the phase-space sequence - is
  // identical to an unpinned run. That is what lets two runs be summed point by
  // point to recover the helicity sum.
  static const int pinned(ATOOLS::Settings::GetMainSettings()["CEEX"]["PIN_PHOTON_HEL"].Get<int>());
  m_PhoHel.clear();
  for (size_t i(0); i < m_allphotons.size(); ++i) {
    const int h(ran->Get() < 0.5 ? 1 : -1);
    m_PhoHel.push_back(pinned ? pinned : h);
  }
}


Ceex_Base::StageTableSwap::StageTableSwap(Ceex_Base &cb) :
  c(cb), legs(cb.m_stagelegs), mom(cb.m_stagemom), res(cb.m_resonances),
  sfac(cb.m_Sfac), stage(cb.m_stage), fixed(cb.m_fixedstage),
  reduces(cb.m_stagereduces), wgroups(cb.m_wgroups), groups(cb.m_groups),
  initstage(cb.m_initstage), nstages(cb.m_nstages) {}

Ceex_Base::StageTableSwap::~StageTableSwap()
{
  c.m_stagelegs = legs; c.m_stagemom = mom; c.m_resonances = res;
  c.m_Sfac = sfac; c.m_stage = stage; c.m_fixedstage = fixed;
  c.m_stagereduces = reduces; c.m_wgroups = wgroups; c.m_groups = groups;
  c.m_initstage = initstage; c.m_nstages = nstages;
}

double Ceex_Base::CrudeFromGenerator()
{
  static const int crudeborn(ATOOLS::Settings::GetMainSettings()["CEEX"]
                             ["CRUDE_BORN"].SetDefault(1).Get<int>());
  static const int pfmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["NO_PSEUDOFLUX"].Get<int>());
  if (crudeborn == 0 || !m_comixborn || !m_cxbalignok
      || m_pceex.size() != m_flavs.size() || m_flavs.size() < 4) return -1.;
  if (m_flavs.size() > 4 && m_prefsr.size() != m_flavs.size()) return -1.;
  if (m_PhoHel.size() != m_allphotons.size()) return -1.;
  StageTableSwap keep(*this);
  // the generator's stages: the flat groups, no resonance
  m_wgroups.clear();
  m_groups = m_crudegroups;
  CalculateSfactors();
  if ((int)m_Sfac.size() != m_nstages
      || (!m_allphotons.empty() && m_Sfac[0].size() != m_allphotons.size()))
    return -1.;
  int last(0);
  PartitionStartCollapsed(last);
  const int nhel(Amplitude::NHel()), fmaskr(Amplitude::NHel() - 1);
  double rho(0.);
  for (;;) {
    Vec4D PX(m_pceex[0] + m_pceex[1]);
    Complex sProd(1., 0.);
    double crudeprod(1.);
    bool crudefixed(false);
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      Complex sj(0., 0.);
      if (m_fixedstage.size() == m_allphotons.size() && m_fixedstage[j]) {
        double inc(0.);
        for (int g(0); g < m_nstages; ++g) {
          sj += m_Sfac[g][j]; inc += std::norm(m_Sfac[g][j]); }
        crudeprod *= inc;
        crudefixed = true;
      } else {
        sj = m_Sfac[m_stage[j]][j];
        crudeprod *= std::norm(sj);
      }
      sProd *= sj;
      if (m_stagereduces[m_stage[j]]) PX -= m_allphotons[j];
    }
    const double X2(PX.Abs2());
    if (X2 > 0. && m_svarQ > 0.) {
      double pflux(1.);
      if (pfmode != 1)
        for (int g(0); g < m_nstages; ++g) {
          if (g == m_initstage || m_stagelegs[g].empty()) continue;
          const Vec4D q(StageSystemMomentum(g)), K(StagePhotonSum(g, -1));
          const double q2(q.Abs2());
          if (q2 > 0.) pflux *= (q + K).Abs2()/q2;
        }
      // the generator's Born at its own point, as InfraredSubtractedME_0_0
      Vec4D_Vector pb;
      Amplitude Cred;
      const bool legs(m_flavs.size() == 4 ? BornLegsAt(PX, pb)
                                           : GeneratorBornAt(PX, pb));
      if (!legs || !ComixBornAmplitude(pb, Cred, NULL, -1., -1.)) return -1.;
      const double crudered((pfmode == 1 ? 1. : pflux) * (m_s/X2)
                            * (crudefixed ? crudeprod : std::norm(sProd)));
      double rc(0.);
      for (int f(0); f < nhel; ++f)
        rc += crudered * std::norm(m_cxbalign.m_A[f]
                                   * Cred.m_A[f ^ (m_comixflip & fmaskr)]);
      rho += rc/4.;
    }
    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
  return rho;
}

/*
  PartitionStart, with the softest photons collapsed onto a fixed stage
  when more than CEEX: MAX_PARTITION_PHOTONS would be enumerated. The
  numerator partition sum and the generator-side crude (CrudeFromGenerator)
  both start from this, on their own stage tables.
*/
void Ceex_Base::PartitionStartCollapsed(int &last) {
  static const size_t maxphot(Settings::GetMainSettings()["CEEX"]["MAX_PARTITION_PHOTONS"]
                              .SetDefault(100).Get<int>());
  if (m_allphotons.size() > maxphot) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): "<<m_allphotons.size()<<" photons would need "
                 <<m_nstages<<"^"<<m_allphotons.size()<<" partitions. Above "
                 <<"CEEX: MAX_PARTITION_PHOTONS = "<<maxphot<<" the softest "
                 <<"photons are collapsed (fixed stage, total eikonal, as "
                 <<"below SOFT_PARTITION_CUT) until that many are left to "
                 <<"enumerate. Reported once."<<std::endl;
    }
    /*
      The earlier fallback put EVERY photon on the initial stage. For a hard
      photon collinear to a final-state lepton that leaves the final-state
      emission part of M_1 unsubtracted against a Born off the resonance:
      measured |beta_1| = 790 |A_0| on a 15-photon event at 250 GeV with a
      28 GeV photon 6 mrad from the muon, a weight of 4e5 times the crude.
      Collapsing the softest photons instead is the same approximation the
      soft cut makes, applied adaptively, and keeps the hard ones enumerated.
    */
    PartitionStart(last);
    size_t nfree(0);
    for (size_t j(0); j < m_allphotons.size(); ++j)
      if (!m_fixedstage[j]) ++nfree;
    while (nfree > maxphot) {
      size_t js(m_allphotons.size()); double emin(-1.);
      for (size_t j(0); j < m_allphotons.size(); ++j)
        if (!m_fixedstage[j] && (emin < 0. || m_allphotons[j][0] < emin))
          { emin = m_allphotons[j][0]; js = j; }
      if (js == m_allphotons.size()) break;
      m_fixedstage[js] = 1;
      int best(m_initstage); double bmax(-1.);
      for (int g(0); g < m_nstages; ++g)
        if (js < m_Sfac[g].size() && std::abs(m_Sfac[g][js]) > bmax)
          { bmax = std::abs(m_Sfac[g][js]); best = g; }
      m_stage[js] = HasFSR() ? best : m_initstage;
      --nfree;
    }
    last = (HasFSR() && nfree > 0) ? 0 : 1;
  } else {
    PartitionStart(last);
  }
}

void Ceex_Base::Calculate() {
  // Lab copies first: everything below boosts into the CEEX frame.
  m_allphotons_lab = m_isrphotons;
  m_allphotons_lab.insert(m_allphotons_lab.end(),
                          m_fsrphotons.begin(), m_fsrphotons.end());
  for (size_t i(0); i < m_isrphotons.size(); ++i) m_cms.Boost(m_isrphotons[i]);
  for (size_t i(0); i < m_fsrphotons.size(); ++i) m_cms.Boost(m_fsrphotons[i]);
  for (size_t i(0); i < m_bornmomenta.size(); ++i) m_cms.Boost(m_bornmomenta[i]);

  m_justdumped = false;
  m_rhocrud = 0.0;
  { static const bool coh(Settings::GetMainSettings()["CEEX"]["IFI"]
                          .SetDefault(1).Get<int>() != 0);
    m_ifi_coherent = coh || m_order == 2; }
  m_inc00 = m_inc11 = m_inc01 = 0.;
  m_b1n = 0; m_b1min = m_b1max = m_b1sum = m_b1sq = 0.;  // beta_1 spread, per event
  m_beta10 = 0.0;
  m_beta01 = 0.0;
  m_beta00 = 0.0;
  BuildCeexMomenta();
  RepairMomentumBalance();
  /*
    Born-configuration probe (env SHERPA_CEEX_BORNNORM), once per event.

  */
  { static const bool bn(ATOOLS::Settings::GetMainSettings()["CEEX"]["BORNNORM_CHECK"].Get<int>() != 0);
    if (bn && m_bornmomenta.size() >= 4 && m_pceex.size() >= 4) {
      const Vec4D Qb(m_bornmomenta[m_if1] + m_bornmomenta[m_if2]);
      const Vec4D Qp(m_pceex[m_if1] + m_pceex[m_if2]);
      Vec4D q2(m_pceex[m_if1]), b0(m_pceex[0]);
      Poincare cm(Qp);
      cm.Boost(q2); cm.Boost(b0);
      const double n1(Vec3D(q2).Abs()), n2(Vec3D(b0).Abs());
      const double cth(n1 > 0. && n2 > 0. ? (Vec3D(q2)*Vec3D(b0))/(n1*n2) : 0.);
      const double sp_save(m_sp);
      m_sp = Qp.Abs2();
      MakeProp();
      Amplitude AB;
      BornAmplitude(m_pceex, AB, -1., -1., -1);
      double sum(0.);
      const int nh(Amplitude::NHel());
      for (int f = 0; f < nh; ++f) {
        const Complex a(m_e * m_e * AB.m_A[f]);
        sum += std::real(a * conj(a));
      }
      m_sp = sp_save;
      MakeProp();
      const double avg(sum / 4.);
      std::cerr << "@@@ BORNANG m_born=" << m_born
                << " sphys=" << Qp.Abs2() << " sborn=" << Qb.Abs2()
                << " cth=" << cth
                << " ceexborn=" << avg
                << " ratio=" << (m_born != 0. ? avg/m_born : 0.) << std::endl;
    } }
  ZerAmplit();
  m_allphotons = m_isrphotons;
  m_allphotons.insert(m_allphotons.end(),
                      m_fsrphotons.begin(), m_fsrphotons.end());
  MakePhotonHel();
  /*
    Frame probe (CEEX: FRAME_PROBE): does what CEEX is handed balance? Beams
    minus final legs minus every photon, and the same with ISR or FSR photons
    only, so a frame mismatch between the legs and one photon family shows
    up as a residual of the size of that family's momentum. Also the angle
    of each FSR photon to its nearest final lepton, to be compared with the
    same angle in the event record.
  */
  { static const int fp(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["FRAME_PROBE"].SetDefault(0).Get<int>());
    static int nfp(0);
    // FRAME_PROBE: 1 prints the first 40 events, N > 1 the first N
    if (fp && nfp < (fp > 1 ? fp : 40) && m_pceex.size() >= 4) {
      ++nfp;
      Vec4D P(m_pceex[0] + m_pceex[1]), KI, KF;
      for (size_t i(2); i < m_pceex.size(); ++i) P -= m_pceex[i];
      for (const Vec4D &k : m_isrphotons) KI += k;
      for (const Vec4D &k : m_fsrphotons) KF += k;
      const Vec4D r(P - KI - KF);
      std::cerr<<std::setprecision(5)<<"@@@ FRAME nI="<<m_isrphotons.size()
               <<" nF="<<m_fsrphotons.size()
               <<" |KI|="<<KI[0]<<" |KF|="<<KF[0]
               <<" resid=("<<r[0]<<","<<r[1]<<","<<r[2]<<","<<r[3]<<")"
               <<" resid_noF3="<<(P-KI)[3]<<" resid_noI3="<<(P-KF)[3]
               <<" Pfin3="<<(m_pceex[0]+m_pceex[1]-P)[3]
               <<" KI3="<<KI[3]<<" KF3="<<KF[3]<<std::endl;
      /*
        How far off its nominal mass shell is each leg Comix will be handed?
        Comix rebuilds every leaf exactly on the flavour's shell (SetPWide)
        and, in the CEEX entry points, does NOT project the set back onto
        momentum conservation; a leg off shell by dm2 moves by dm2/(2p) and
        the complement-route propagators inherit that as an ABSOLUTE
        virtuality error, fatal for a photon inside sqrt(E/E_k) m/E of a
        beam.
      */
      std::cerr<<std::setprecision(3)<<"@@@ FRAME shell";
      for (size_t i(0); i < m_pceex.size(); ++i)
        std::cerr<<" ["<<i<<"] "<<(m_pceex[i].Abs2() - sqr(m_flavs[i].Mass()));
      for (const Vec4D &k : m_isrphotons) std::cerr<<" I "<<k.Abs2();
      for (const Vec4D &k : m_fsrphotons) std::cerr<<" F "<<k.Abs2();
      std::cerr<<" |resid|="<<std::setprecision(3)
               <<sqrt(sqr(r[0])+sqr(r[1])+sqr(r[2])+sqr(r[3]))
               <<" before_repair="<<m_repairresid<<std::endl;
      /*
        At one photon: the generator's Born point (m_bornmomenta) against the
        reduced legs CEEX builds for the same partition (BornLegsAt(P - k)).
        The crude Born the fixed-order weight divides by sits at the former;
        the partition crude at the latter.
      */
      if (m_isrphotons.size() + m_fsrphotons.size() == 1 && m_bornmomenta.size() >= 4) {
        Vec4D_Vector pb;
        const Vec4D k(m_isrphotons.empty() ? m_fsrphotons[0] : m_isrphotons[0]);
        if (BornLegsAt(m_pceex[0] + m_pceex[1] - k, pb)) {
          std::cerr<<std::setprecision(6)<<"@@@ FRAME born1";
          for (size_t i(0); i < 4 && i < pb.size(); ++i)
            std::cerr<<" gen["<<i<<"]="<<m_bornmomenta[i]<<" red["<<i<<"]="<<pb[i];
          std::cerr<<std::endl;
          // the Comix Born squared at both points, against the fixed-order
          // crude Born m_born the handler divides by
          Amplitude Ag, Ar;
          Vec4D_Vector pg(m_bornmomenta.begin(), m_bornmomenta.begin() + Min(m_bornmomenta.size(), m_pceex.size()));
          double ng(0.), nr(0.);
          if (pg.size() == m_pceex.size() && ComixBornAmplitude(pg, Ag, NULL, -1., -1.))
            for (int f(0); f < Amplitude::NHel(); ++f) ng += std::norm(Ag.m_A[f]);
          if (ComixBornAmplitude(pb, Ar, NULL, -1., -1.))
            for (int f(0); f < Amplitude::NHel(); ++f) nr += std::norm(Ar.m_A[f]);
          const double x(2.*k[0]/sqrt(m_s));
          // and at the generator's own s' point, m_plabmom (lab frame; the
          // Born is Lorentz invariant, so no boost is needed for the norm)
          double np(-1.);
          if (m_plabmom.size() == m_pceex.size()) {
            Amplitude Ap;
            if (ComixBornAmplitude(m_plabmom, Ap, NULL, -1., -1.)) {
              np = 0.; for (int f(0); f < Amplitude::NHel(); ++f) np += std::norm(Ap.m_A[f]); }
            std::cerr<<"@@@ FRAME plab";
            for (size_t i(0); i < m_plabmom.size(); ++i) std::cerr<<" ["<<i<<"]="<<m_plabmom[i];
            std::cerr<<std::endl;
          }
          std::cerr<<"@@@ FRAME bornsq x="<<x<<" m_born="<<m_born
                   <<" B2plab/m_born="<<(m_born!=0.?np/m_born:-1.)
                   <<" B2gen="<<ng<<" B2red="<<nr
                   <<" B2gen/m_born="<<(m_born!=0.?ng/m_born:-1.)
                   <<" B2red/m_born="<<(m_born!=0.?nr/m_born:-1.)
                   <<" B2red/m_born*(1-x)="<<(m_born!=0.?nr/m_born*(1.-x):-1.)<<std::endl;
        }
      }
    } }
  CalculateSfactors();
  /*
    CalculateSfactors clears m_Sfac when BuildStages rejects the decomposition
    (a stage carrying net charge, so gauge dependent). Bail here rather than
    index an empty table below: the event then contributes no CEEX weight,
    which YFS_Handler already treats as "CEEX produced nothing".
  */
  if ((int)m_Sfac.size() != m_nstages
      || (!m_allphotons.empty() && m_Sfac[0].size() != m_allphotons.size())) {
    m_nparts = 0;
    return;
  }
  m_spincache.resize(1 + 4*m_allphotons.size());
  m_spinvalid.assign(m_spincache.size(), 0);
  m_realphot.assign(m_allphotons.size(), Amplitude());
  m_realphotM1.assign(RealVirtualMode() == 2 ? m_allphotons.size() : 0, Amplitude());
  /*
    Trace the beta_1 pieces on the event the CHECK_XS dump will write - the
    same gate the dump uses, evaluated before the partition loop so the loop
    can print as it goes. One event, so a 2^n fan-out is affordable.
  */
  { static const bool tr(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["BETA1_TRACE"].Get<int>() != 0);
    static const double dx(ATOOLS::Settings::GetMainSettings()["CEEX"]
                           ["DUMP_XMIN"].Get<double>());
    static const size_t dn(ATOOLS::Settings::GetMainSettings()["CEEX"]
                           ["DUMP_NPHOT"].Get<int>());
    const double xg(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                    2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass()
                    : 0.);
    static const double heavy(ATOOLS::Settings::GetMainSettings()["CEEX"]
                              ["TRACE_FACTOR_ABOVE"].Get<double>());
    m_b1trace = (tr && m_checkxs && !m_ceexdumped
                 && m_allphotons.size() == dn && xg > dx)
                || (heavy > 0. && !m_allphotons.empty()); }
  /*
    The invariant of the outgoing FERMION pair, not of legs 2 and 3. With
    anything else in the final state those are not the same: for
    e+e- -> H l+ l- leg 2 is the Higgs, so this was (p_H + p_l-)^2 rather than
    the dilepton mass. m_svarQ is the denominator of the pseudo-flux factor
    m_cfac = sProd*(svarX/m_svarQ) and feeds the FSR kinematic term, so it is
    load bearing, and being wrong there is silent - it is a perfectly good
    invariant, just not the one the decay propagator sits at.
  */
  m_svarQ = (m_pceex[m_if1] + m_pceex[m_if2]).Abs2();
  m_sQ    = m_svarQ;
  MakePropT(m_pceex);

  
  // One Comix/hand ratio per photon, before the partition loop: inside it the
  // soft weights are in scope but the photon's TOTAL amplitude is not.
  { static const bool bs(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["BETA1_SCAN"].Get<int>() != 0);
    if (bs) Beta1Scan(); }
  BuildComixBornAlignment();
  BuildComixPhotonRatios();
  // CEEX: TCHANNEL_REDUCED_BORN (-1 auto, see Ceex_Base::RegisterDefaults)
  { static const int rb(ATOOLS::Settings::GetMainSettings()["CEEX"]
                        ["TCHANNEL_REDUCED_BORN"].Get<int>());
    /*
      Auto (-1) needs a single radiating stage as well as an exchange line.
      ComixBeta1At subtracts ONE reduced Born, times the soft factors of all
      stages, while the partition sum puts a different reduced Born (at that
      partition's X) on each stage, so rho_1 = |M_1|^2 at one photon holds
      only when there is one stage. e+e- -> gamma gamma has one; Bhabha has a
      final-state stage too. Measured on Bhabha at the Z pole, Born+real,
      2026-09-26: with 1, CEEX/YFS.NLO = 0.992 (one ISR photon), 1.055 (one
      FSR photon), 1.145 (two or more), column +10.6%; with 0, 0.992, 0.997,
      0.991, column +0.3%.
    */
    // Radiating stages only: a neutral final state still gets an EMPTY
    // placeholder stage (BuildStages), so gamma gamma has m_nstages = 2.
    int nrad(0);
    for (size_t g(0); g < m_stagelegs.size(); ++g)
      if (!m_stagelegs[g].empty()) ++nrad;
    m_redborn = rb > 0 || (rb < 0 && m_comixborn && nrad == 1
                           && BornHasExchangeLine());
    // the reduced point is LegsAt's 2 -> 2 construction; with W stages the
    // partition Born stays the shifted physical-spinor one
    if (WStagesActive() && m_redborn) {
      static bool warned(false);
      if (!warned) { warned = true;
        msg_Error()<<METHOD<<"(): CEEX: TCHANNEL_REDUCED_BORN is not defined "
                   <<"with W stages; ignored."<<std::endl; }
      m_redborn = false;
    } }

  int last(0), nparts(0);
  PartitionStartCollapsed(last);
  if (m_allphotons.size() > m_maxnphot) m_maxnphot = m_allphotons.size();
  // CEEX: ORDER 2 - which photons can form beta_2 pairs (after the collapse)
  const bool dobeta2(m_order == 2 && PrepareBeta2());
  if (dobeta2) ++m_b2events;
  Amplitude snap0, snap1;
  for (;;) {
    ++nparts;
    if (!m_ifi_coherent) { snap0 = m_AmpExpo0; snap1 = m_AmpExpo1; }
    Vec4D PX(m_pceex[0] + m_pceex[1]);
    Complex sProd(1., 0.);
    double crudeprod(1.);
    bool crudefixed(false);
    std::vector<Complex> Sactu(m_allphotons.size(), Complex(1., 0.));
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      // an enumerated photon carries its stage's eikonal, a fixed (soft)
      // one the total: its sum over stages has been done
      if (m_fixedstage.size() == m_allphotons.size() && m_fixedstage[j]) {
        Sactu[j] = Complex(0., 0.);
        double inc(0.);
        for (int g(0); g < m_nstages; ++g) {
          Sactu[j] += m_Sfac[g][j];
          inc += std::norm(m_Sfac[g][j]);
        }
        // the crude is the generator's density: incoherent over the labels
        crudeprod *= inc;
        crudefixed = true;
      } else {
        Sactu[j] = m_Sfac[m_stage[j]][j];
        crudeprod *= std::norm(Sactu[j]);
      }
      sProd   *= Sactu[j];
      if (m_stagereduces[m_stage[j]]) PX -= m_allphotons[j];
    }
    m_sactu = Sactu;
    m_crudeprod = crudeprod;
    m_crudefixed = crudefixed;
    /*
      W_EIKONAL 1: the resonance stages' soft factors depend on THIS
      partition (the W's pole momentum carries the partition's decay
      photons), so they are rebuilt here and the products above redone
      with them. No-op without W stages.
    */
    if (WStagesActive() && m_weikonal != 0) {
      RecomputeResonanceSfactors();
      sProd = Complex(1., 0.); crudeprod = 1.;
      for (size_t j(0); j < m_allphotons.size(); ++j) {
        if (m_fixedstage.size() == m_allphotons.size() && m_fixedstage[j]) {
          Sactu[j] = Complex(0., 0.);
          double inc(0.);
          for (int g(0); g < m_nstages; ++g) {
            Sactu[j] += m_Sfac[g][j]; inc += std::norm(m_Sfac[g][j]); }
          crudeprod *= inc;
        } else {
          Sactu[j] = m_Sfac[m_stage[j]][j];
          crudeprod *= std::norm(Sactu[j]);
        }
        sProd *= Sactu[j];
      }
      m_sactu = Sactu;
      m_crudeprod = crudeprod;
    }
    const double svarX(PX.Abs2());
    if (m_checkxs) {
      const Complex pz(1./Complex(svarX - sqr(m_MZ), m_gZ*svarX/m_MZ));
      const double a(std::abs(pz));
      if (nparts == 1) { m_pzmin = m_pzmax = a; }
      else { m_pzmin = Min(m_pzmin, a); m_pzmax = Max(m_pzmax, a); }
    }
    if (svarX <= 0. || m_svarQ <= 0.) {
      if (last) break;
      PartitionPlus(last);
      if (last == 2) break;
      continue;
    }
    m_sp   = svarX;
    m_Sprod = sProd;
    /*
      The pseudo-flux factor relates the invariant flowing INTO the radiating
      final-state system to that system's own invariant, m_svarQ. KKMC writes
      it as svarX/svarQ (KKceex.cxx:395) and hep-ph/0006359 eq.(42) fixes the
      meaning: it "disappears in the in-space situation p_a+p_b = p_c+p_d" and
      "really matters if at least one hard FSR photon is present" - i.e. it is
      1 when there is no FSR photon.

      svarX is the whole X system. At 2 -> 2 that IS the radiating system
      (X = q_c + q_d + sum_FSR k by momentum conservation), so svarX/svarQ is
      1 + O(x_gamma) and cancels against the FSR (1-CKine) term below, which is
      what KKMC designed it to do ("Contribution -svarX/svarQ from HERE cancels
      exactly with svarX/svarQ in beta0", KKceex.cxx:1499).

      Past 2 -> 2 the two part company: for e+e- -> H mu+ mu- the numerator
      carries the Higgs and the denominator does not, so the ratio is s/m_ll^2
      ~ 7.5 with NO photon at all. Squared that is 56.5, which is exactly the
      factor by which the Comix real was observed to suppress the CEEX weight -
      the real carries no m_cfac while m_rhocrud carries m_cfac^2.

      So subtract the final legs that are NOT the radiating pair. At 2 -> 2 the
      loop body never runs, PY is PX itself and svarY is bit-for-bit svarX, so
      the 2 -> 2 stream is untouched by construction rather than by rounding.
    */
    Vec4D PY(PX);
    for (size_t i(2); i < m_pceex.size(); ++i)
      if (i != m_if1 && i != m_if2) PY -= m_pceex[i];
    const double svarY(PY.Abs2());
    m_svarY = svarY;          // the decay line's own partition-shifted scale
    m_PXvec = PX;             // and the momentum that realises m_sp
    /*
      The pseudo-flux is KKMC bookkeeping, not physics on its own: it exists so
      that svarX/svarQ in beta_0 cancels the (1-CKine) term in the hand-coded
      FSR beta_1 (KKceex.cxx:1499, and the comment above). A beta_1 taken
      straight from the amplitude has no CKine term for it to cancel against,
      so the factor is left dangling on beta_0 alone. NO_PSEUDOFLUX drops it
      from both sides at once, which is the only self-consistent way to have it
      off.
    */
    /*
      NO_PSEUDOFLUX: 0 = svarY/svarQ on beta_0 in both rho_0 and rho_1 (KKMC's
      beta_0); 1 = on neither; 2 = on rho_0 only. KKMC's O(alpha^1) AMPLITUDE
      is flux-free up to 2k_i.k_j/Q^2 - it adds sProd(1-CKine)B per
      final-state photon (KKceex.cxx:1552,1569) against the flux on beta_0 -
      while its RhoExp0 keeps the flux, so 2 is KKMC's own rho_1/rho_0.
      Measured on the seed-11 n=2 point: rho_1 with the flux off agrees with
      KKMC's to 0.04%; with it on, A1 loses 0.26 A0 in every helicity.
    */
    static const int pfmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                            ["NO_PSEUDOFLUX"].Get<int>());
    /*
      The pseudo-flux, generalised: one factor (q_g + K_g)^2 / q_g^2 per
      final stage g, q_g the physical momentum of the stage's legs and K_g
      the photons this partition assigns to it. At 2 -> 2 q + K_F = P - K_I
      by conservation, so this is KKMC's svarX/svarQ to rounding; for H l l
      it is the (X - p_H)^2 / m_ll^2 the code carried by hand; for two
      resonances each gets its own. It stays in rho_0 only (mode 2): the
      2 -> 2 comparison with KKMC fixes that, and for H mu mu dropping it
      moved CEEX from -6.4% to +11.3% against Born+real.
    */
    m_pflux = 1.;
    if (pfmode != 1)
      for (int g(0); g < m_nstages; ++g) {
        if (g == m_initstage || m_stagelegs[g].empty()) continue;
        // for a W decay stage q is the W's daughters: (P_W + K)^2/P_W^2
        const Vec4D q(StageSystemMomentum(g)), K(StagePhotonSum(g, -1));
        const double q2(q.Abs2());
        if (q2 > 0.) m_pflux *= (q + K).Abs2()/q2;
      }
    m_cfac = sProd * Complex(pfmode == 0 ? m_pflux : 1., 0.);   // rho_1 side
    MakeProp();
    // Electroweak form factors at THIS partition's scale, and the scattering
    // angle the WW/ZZ boxes depend on. No-op unless CEEX: WEAK is set.
    MakeEWFF(m_sp, cos(m_pceex[m_if1].Theta()));
    MakeBoxMandelstams(PX);

    /*
      Does Comix's beta_1 depend on the partition? If it does not, one
      evaluation pair per photon serves all 2^n partitions and the per-photon
      real is affordable; if it does, the cost is 2 x n x 2^n per event.
    */
    { static const bool b1chk(ATOOLS::Settings::GetMainSettings()["CEEX"]
                              ["BETA1_PARTITION_CHECK"].Get<int>() != 0);
      if (b1chk && !m_allphotons.empty()) {
        Amplitude B1;
        if (ComixBeta1At(m_allphotons[0], m_PhoHel[0], m_sp, B1)) {
          double n2(0.);
          for (int f(0); f < Amplitude::NHel(); ++f) n2 += std::norm(B1.m_A[f]);
          const double nb(sqrt(n2));
          if (nb > 0.) {
            if (m_b1n == 0) { m_b1min = m_b1max = nb; }
            else { m_b1min = Min(m_b1min, nb); m_b1max = Max(m_b1max, nb); }
            m_b1sum += nb; m_b1sq += nb*nb;
            ++m_b1n;
          }
        }
      } }
    if (dobeta2) ResetBeta1Record();
    InfraredSubtractedME_0_0();
    m_b1terms.clear();        // CEEX: TCHANNEL_MULTIPHOTON 2: this partition's
    // Only when a virtual was asked for and CEEX is its source; see
    // Ceex_Base::CeexOwnVirtual.
    if (CeexOwnVirtual()) InfraredSubtractedME_0_1();
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      if (m_checkxs && m_stagereduces[m_stage[j]]) {
        CheckSoftFactor(m_allphotons[j]);
        CheckUVDiagonality(m_allphotons[j]);
      }
      /*
        beta_1: Comix's one-photon amplitude with the eikonal subtracted off
        it. There is no hand-coded fallback. The KKMC-shaped routines it
        replaced are a formula for e+e- -> f fbar, not a method - they put the
        photon's spinor index into a beam slot, which only means anything for a
        single s-channel 2 -> 2 - and for e+e- -> nu nu they return a CEEX
        cross section 3.1x the NLO. Falling back to them hides a wrong answer
        behind a correct-looking one, so a failure here contributes nothing and
        says so.
      */
      if (!ComixInfraredSubtracted_1_0(m_allphotons[j], m_PhoHel[j],
                                       (int)j)) {
        static bool warned(false);
        if (!warned) {
          warned = true;
          msg_Error()<<METHOD<<"(): beta_1 unavailable from Comix; it is "
                     <<"being dropped, not replaced.\n";
        }
      }
    }
    // CEEX: TCHANNEL_MULTIPHOTON 2: the factorised remainder of the beta_1
    // terms just added (empty, so a no-op, in every other mode)
    AddFactorisedRemainder();
    /*
      CEEX: ORDER 2 - beta_2 for every eligible pair, after the beta_1 loop
      of the SAME partition, whose M_1 and eikonals it subtracts.
    */
    if (dobeta2)
      for (size_t j(0); j < m_allphotons.size(); ++j)
        for (size_t l(j + 1); l < m_allphotons.size(); ++l)
          if (Beta2PairEligible((int)j, (int)l))
            ComixInfraredSubtracted_2_0((int)j, (int)l);
    { static const int kkflux(ATOOLS::Settings::GetMainSettings()["CEEX"]
                              ["KKMC_FLUX_EMULATION"].Get<int>());
      if (kkflux && m_order == 2 && m_flavs.size() == 4) {
        const Vec4D Q(m_pceex[m_if1] + m_pceex[m_if2]);
        const double Q2(Q.Abs2());
        double rall(0.), rsoft(0.);
        for (size_t j(0); j < m_allphotons.size(); ++j)
          for (size_t l(j + 1); l < m_allphotons.size(); ++l)
            if (!m_stagereduces[m_stage[j]] && !m_stagereduces[m_stage[l]]) {
              const double t(2.*(m_allphotons[j]*m_allphotons[l])/Q2);
              rall += t;
              if (!Beta2PairEligible((int)j, (int)l)) rsoft += t;
            }
        for (int f(0); f < Amplitude::NHel(); ++f) {
          m_AmpFluxAll.m_A[f]  += m_Sprod*rall*m_partborn0.m_A[f];
          m_AmpFluxSoft.m_A[f] += m_Sprod*rsoft*m_partborn0.m_A[f];
        }
      } }

    if (!m_ifi_coherent)   // this partition's own contribution, squared
      for (int f(0); f < Amplitude::NHel(); ++f) {
        const Complex d0(m_AmpExpo0.m_A[f] - snap0.m_A[f]);
        const Complex d1(m_AmpExpo1.m_A[f] - snap1.m_A[f]);
        m_inc00 += std::norm(d0); m_inc11 += std::norm(d1);
        m_inc01 += std::real(std::conj(d0)*d1);
      }
    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
  /*
    CEEX: CRUDE_FROM_GENERATOR - the crude on the generator's own stages.
    Mode 2 only prints it next to the per-partition one (the gate: the two
    must agree wherever CEEX's stages are the generator's dipoles).
  */
  const int crudegen(m_crudegen >= 0 ? m_crudegen : (WStagesActive() ? 1 : 0));
  if (crudegen != 0) {
    const double rg(CrudeFromGenerator());
    if (crudegen == 2) {
      static long ncg(0);
      if (ncg < 5000) { ++ncg;
        std::cerr<<std::setprecision(10)<<"@@@ CRUDEGEN nphot="<<m_allphotons.size()
                 <<" nstages="<<m_nstages<<" rhocrud="<<m_rhocrud<<" gen="<<rg
                 <<" ratio="<<(m_rhocrud > 0. && rg >= 0. ? rg/m_rhocrud : -1.)
                 <<std::endl; }
    } else if (rg >= 0.) m_rhocrud = rg;
    else {
      static bool warned(false);
      if (!warned) { warned = true;
        msg_Error()<<METHOD<<"(): CEEX: CRUDE_FROM_GENERATOR asked for but the "
                   <<"generator crude could not be built; the per-partition "
                   <<"crude is used. Reported once."<<std::endl; }
    }
  }
  if (m_b1trace) {
    double n0(0.), n1(0.);
    for (int f(0); f < Amplitude::NHel(); ++f) {
      n0 += std::norm(m_AmpExpo0.m_A[f]); n1 += std::norm(m_AmpExpo1.m_A[f]); }
    // per helicity, in the same digit order as the KKMC harness's KKPARTP
    for (int j1 = 0; j1 <= 1; ++j1)
      for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
          for (int j4 = 0; j4 <= 1; ++j4) {
            const Complex a0(m_AmpExpo0.m_A[Idx(j1,j2,j3,j4)]);
            const Complex a1(m_AmpExpo1.m_A[Idx(j1,j2,j3,j4)]);
            if (std::abs(a0) < 1e-6*sqrt(n0) && std::abs(a1) < 1e-6*sqrt(n1)) continue;
            std::cerr<<std::setprecision(8)<<"@@@ B1HEL "<<j1<<j2<<j3<<j4
                     <<" A0="<<a0<<" A1-A0="<<(a1-a0)
                     <<" |A1-A0|/|A0|="<<(std::abs(a0)>0.? std::abs(a1-a0)/std::abs(a0) : -1.)
                     <<" arg((A1-A0)/A0)="<<(std::abs(a0)>0.? std::arg((a1-a0)/a0) : 0.)
                     <<std::endl;
          }
    std::cerr<<"@@@ B1SUM |A0|="<<sqrt(n0)<<" |A1|="<<sqrt(n1)
             <<" nparts="<<nparts;
    for (size_t j(0); j < m_realphot.size(); ++j) {
      double nb(0.);
      for (int f(0); f < Amplitude::NHel(); ++f) nb += std::norm(m_realphot[j].m_A[f]);
      std::cerr<<" |beta1_"<<j<<"|="<<sqrt(nb)
               <<" E_"<<j<<"="<<m_allphotons[j][0];
    }
    std::cerr<<std::endl;
  }
  /*
    Closure at one photon. The partition sum ought to return Comix's exact
    one-photon amplitude, because the stage weights sum to one. It cannot be
    exact: the I and F terms carry the propagator at X_I and X_F, and that
    difference IS the split between photons off the initial and final legs -
    the eikonal weights only approximate it. So this number is the size of
    that approximation, measured against the exact amplitude rather than
    against another reconstruction.
  */
  { static const bool cl(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["CLOSURE_1PHOT"].Get<int>() != 0);
    if (cl && m_allphotons.size() == 1 && m_cxbalignok) {
      Amplitude Mex;
      Vec4D_Vector pp(m_pceex);
      pp.push_back(m_allphotons[0]);
      const double rn(RealNorm());
      if (rn > 0. && ComixRealAt(pp, m_PhoHel[0], Mex, -1.)) {
        double nd(0.), ne(0.);
        for (int f(0); f < Amplitude::NHel(); ++f) {
          const Complex ex(m_cxbalign.m_A[f]*Mex.m_A[f]/rn);
          nd += std::norm(m_AmpExpo1.m_A[f] - ex);
          ne += std::norm(ex);
        }
        if (ne > 0.) {
          double nb(0.), nr(0.), n1(0.), nbr(0.);
          for (int f(0); f < Amplitude::NHel(); ++f) {
            nb += std::norm(m_snapBorn.m_A[f]); nr += std::norm(m_snapReal.m_A[f]);
            n1 += std::norm(m_AmpExpo1.m_A[f]);
            nbr += std::norm(m_snapBorn.m_A[f] + m_snapReal.m_A[f]); }
          std::cerr<<"@@@ CLOS1 xg="
                   <<(m_s>0.?2.*m_allphotons[0][0]/sqrt(m_s):-1.)
                   <<" rel="<<sqrt(nd/ne)
                   <<" |Born|2="<<nb<<" |Real|2="<<nr<<" |Born+Real|2="<<nbr
                   <<" |A1|2="<<n1<<" |Mex|2="<<ne<<" nparts="<<m_nparts<<std::endl; }
      }
    } }
  { static const bool b1chk(ATOOLS::Settings::GetMainSettings()["CEEX"]
                            ["BETA1_PARTITION_CHECK"].Get<int>() != 0);
    if (b1chk && m_b1n > 1 && m_b1min > 0.)
      { const double mean(m_b1sum/m_b1n);
        const double var(m_b1sq/m_b1n - mean*mean);
        std::cerr<<"@@@ B1PART nphot="<<m_allphotons.size()
                 <<" nparts="<<m_b1n
                 <<" spread="<<(m_b1max/m_b1min)
                 <<" maxovmean="<<(mean>0.?m_b1max/mean:-1.)
                 <<" cv="<<(mean>0.?sqrt(var>0.?var:0.)/mean:-1.)
                 <<std::endl; } }
  /*
    Comix's beta_1 against what the partition loop ACTUALLY accumulated for
    each photon. m_realphot[j] is that object - the stage pieces with their
    partition weights - whereas HandOnePhotonAmplitude is a standalone
    reconstruction, and a reconstruction is not the reference: it produced
    |beta_1| LARGER than the full amplitude on some photons, which cannot be.
  */
  { static const bool b1(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["BETA1_CHECK"].Get<int>() != 0);
    static int nd(0);
    if (b1 && nd < 14)
      for (size_t j(0); j < m_allphotons.size() && nd < 14; ++j) {
        if (j >= m_cxratiook.size() || !m_cxratiook[j]) continue;
        ++nd;
        double nr(0.), ncs(0.);
        const int nh(Amplitude::NHel());
        for (int f = 0; f < nh; ++f) {
          nr  += std::norm(m_realphot[j].m_A[f]);
          ncs += std::norm(m_cxcsub[j].m_A[f]);
        }
        std::cerr<<"@@@ BETA1REF Ek="<<m_allphotons[j][0]
                 <<" |realphot|="<<sqrt(nr)
                 <<" |Csub|="<<sqrt(ncs)
                 <<" ratio="<<(nr>0.? sqrt(ncs/nr) : -1.)<<std::endl;
      } }
  /*
    Report the per-event spread of V(h) across partitions, then reset. Only
    the helicities that carry weight are shown; the mass-suppressed ones have
    a Born near zero and their ratio means nothing.
  */
  { static const bool vp(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["VIRT_PARTITION_CHECK"].Get<int>() != 0);
    static int nv(0);
    if (vp && nv < 8 && !m_vpmin.empty()) {
      ++nv;
      double worst(0.); int nl(0);
      for (size_t f(0); f < m_vpmin.size(); ++f) {
        if (m_vpmax[f] < m_vpmin[f]) continue;
        ++nl;
        const double mid(0.5*(m_vpmax[f] + m_vpmin[f]));
        if (mid > 0.) worst = Max(worst, (m_vpmax[f] - m_vpmin[f])/mid);
      }
      // the actual range, not just the relative spread: a spread of 2 means
      // the minimum reached zero, which needs to be visible
      double gmin(1e30), gmax(-1e30);
      for (size_t f(0); f < m_vpmin.size(); ++f) {
        if (m_vpmax[f] < m_vpmin[f]) continue;
        gmin = Min(gmin, m_vpmin[f]); gmax = Max(gmax, m_vpmax[f]);
      }
      std::cerr<<"@@@ VIRTPART nparts="<<nparts<<" nphot="<<m_allphotons.size()
               <<" nlive="<<nl<<" worst_spread="<<worst
               <<" |V| in ["<<gmin<<","<<gmax<<"]"
               <<" sQ="<<m_svarQ<<std::endl;
    }
    m_vpmin.clear(); m_vpmax.clear(); }
  m_nparts = nparts;
  if (m_checkxs && nparts > 1 && m_pzmin > 0.) {
    const size_t n(Min(m_allphotons.size(), size_t(15)));
    m_pzspread[n] = Max(m_pzspread[n], m_pzmax/m_pzmin);
    m_pzsum[n] += m_pzmax/m_pzmin;
    ++m_pzn[n];
    const double ecm((m_momenta[0]+m_momenta[1]).Mass());
    for (int c(0); c < 4; ++c) {
      const double xcut(m_pzxcut[c]);
      Vec4D Xhard(m_pceex[0] + m_pceex[1]), Xall(Xhard);
      size_t nsoft(0);
      for (size_t j(0); j < m_allphotons.size(); ++j) {
        const double x(2.*m_allphotons[j][0]/ecm);
        if (x > xcut) { Xhard -= m_allphotons[j]; Xall -= m_allphotons[j]; }
        else { Xall -= m_allphotons[j]; ++nsoft; }
      }
      if (nsoft == 0) continue;
      const double sh(Xhard.Abs2()), sa(Xall.Abs2());
      if (sh <= 0. || sa <= 0.) continue;
      const double a(std::abs(1./Complex(sh - sqr(m_MZ), m_gZ*sh/m_MZ)));
      const double b(std::abs(1./Complex(sa - sqr(m_MZ), m_gZ*sa/m_MZ)));
      const double r(a > b ? a/b : b/a);
      m_pzsoft[c] = Max(m_pzsoft[c], r);
      m_pzsoftsum[c] += r; ++m_pzsoftn[c];
      m_pzsoftsav[c] += double(1u << Min(nsoft, size_t(20)));
    }
  }
 
  /*
    Partition-count sanity check: m_nstages^(photons ENUMERATED). Photons
    below SOFT_PARTITION_CUT, and the softest ones above the
    MAX_PARTITION_PHOTONS cap, sit on a fixed stage (m_fixedstage) and are
    not enumerated, so the count is not 2^n. The earlier form of this check
    predated both the collapse and general stages and fired on every run
    with a soft photon (e.g. "2 partitions for 3 photons, expected 8").
  */
  {
    const size_t n(m_allphotons.size());
    size_t nfree(0);
    if (m_fixedstage.size() == n)
      for (size_t j(0); j < n; ++j) if (!m_fixedstage[j]) ++nfree;
    long want(1);
    if (HasFSR() && m_nstages > 1) {
      for (size_t j(0); j < nfree && want <= (1L << 40); ++j) want *= m_nstages;
      if (want > (1L << 40)) want = nparts;       // overflow guard only
    }
    if (nparts != want) {
      ++m_partbad;
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): partition counter enumerated "<<nparts
                   <<" partitions for "<<n<<" photons ("<<nfree
                   <<" enumerated over "<<m_nstages<<" stages), expected "
                   <<want<<". Reported once."<<std::endl;
      }
    }
    ++m_partn;
  }

  // Comix-supplied O(alpha) real, if CEEX: COMIX_REAL asked for it. It
  // OVERWRITES m_AmpExpo1 for one-photon events, so it has to run after the
  // partition loop has finished accumulating and before MakeRho squares.
  ApplyComixReal();

  MakeRho();
  if (m_order == 2) {
    static const bool b2c(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["BETA2_CLOSURE"].Get<int>() != 0);
    static const bool b2s(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["BETA2_SOFT_TEST"].Get<int>() != 0);
    if (b2c) Beta2Closure();
    if (b2s) Beta2SoftLimit();
  }
  { static const bool bc(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["BETA1_CLOSURE"].Get<int>() != 0);
    if (bc) Beta1Closure(); }
  { static const bool sl(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["SOFT_LIMIT_TEST"].Get<int>() != 0);
    if (sl) SoftLimitTest(); }
  { static const bool sn(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["SOFT_NORM_CHECK"].Get<int>() != 0);
    if (sn && !m_allphotons.empty()) {
      size_t js(0);
      for (size_t i(1); i < m_allphotons.size(); ++i)
        if (m_allphotons[i][0] < m_allphotons[js][0]) js = i;
      SoftNormCheck(m_allphotons[js]); } }

  static const double dumpxmin(ATOOLS::Settings::GetMainSettings()["CEEX"]["DUMP_XMIN"].Get<double>());
  static const size_t dumpn(ATOOLS::Settings::GetMainSettings()["CEEX"]["DUMP_NPHOT"].Get<int>());
  const double xgam(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                    2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass() : 0.);
  // DUMP_XMIN_EACH: every photon above this x (a multi-hard-photon point)
  static const double dumpxeach(ATOOLS::Settings::GetMainSettings()["CEEX"]["DUMP_XMIN_EACH"].Get<double>());
  bool dumpeach(true);
  if (dumpxeach > 0. && m_momenta.size() >= 2)
    for (size_t i(0); i < m_allphotons.size(); ++i)
      if (!(2.*m_allphotons[i][0]/(m_momenta[0]+m_momenta[1]).Mass() > dumpxeach))
        dumpeach = false;
  if (m_checkxs && !m_ceexdumped && m_allphotons.size() == dumpn && dumpeach &&
      xgam > dumpxmin && m_momenta.size() >= 6 && m_result0 != 0.) {
    
    Vec4D bal(m_momenta[0] + m_momenta[1] - m_momenta[4] - m_momenta[5]);
    for (size_t i(0); i < m_allphotons.size(); ++i) bal -= m_allphotons[i];
    const double worst(Max(Max(dabs(bal[0]), dabs(bal[1])),
                           Max(dabs(bal[2]), dabs(bal[3]))));
    if (worst > 1e-6) {
      msg_Error()<<METHOD<<"(): CEEX dump skipped, momentum imbalance "
                 <<worst<<std::endl;
      return;
    }
    m_ceexdumped = true;
    m_justdumped = true;
    std::ofstream pt("ceex_point.dat");
    pt << std::setprecision(17);
    const int idx[5] = {0, 1, 4, 5, -1};
    for (int n(0); n < 4; ++n)
      pt << m_momenta[idx[n]][0] << " " << m_momenta[idx[n]][1] << " "
         << m_momenta[idx[n]][2] << " " << m_momenta[idx[n]][3] << "\n";
    
    pt << m_allphotons[0][0] << " " << m_allphotons[0][1] << " "
       << m_allphotons[0][2] << " " << m_allphotons[0][3] << "\n";
    pt << "NPHOT " << m_allphotons.size() << "\n";
    for (size_t i(0); i < m_allphotons.size(); ++i)
      pt << m_allphotons[i][0] << " " << m_allphotons[i][1] << " "
         << m_allphotons[i][2] << " " << m_allphotons[i][3] << " "
         << m_PhoHel[i] << "\n";
    
    {
      Amplitude bref;
      BornAmplitude(m_pceex, bref);
      msg_Error() << std::setprecision(10);
      msg_Error() << "@@@ SHMASS nflav=" << m_flavs.size()
                << " m1=" << (m_flavs.size()>0? m_flavs[0].Mass():-1.)
                << " m2=" << (m_flavs.size()>1? m_flavs[1].Mass():-1.)
                << " m3=" << (m_flavs.size()>2? m_flavs[2].Mass():-1.)
                << " m4=" << (m_flavs.size()>3? m_flavs[3].Mass():-1.)
                << " mass_I=" << m_mass_I << " mass_F=" << m_mass_F
                << std::endl;
      
      msg_Error() << std::setprecision(16) << "@@@ SHPMASS";
      for (size_t i = 0; i < m_pceex.size(); ++i)
        msg_Error() << " p" << i << "=" << m_pceex[i].Mass();
      msg_Error() << std::setprecision(10) << std::endl;
      msg_Error() << "@@@ SHCPL sin2tw=" << m_sin2tw << " norm=" << m_norm
                << " ve=" << m_ve << " vf=" << m_vf
                << " ae=" << m_ae << " af=" << m_af
                << " qe=" << m_qe << " qf=" << m_qf
                << " CZ1=" << CouplingZ(1., 1) << " CZ0=" << CouplingZ(1., 0)
                << " CG=" << CouplingG() << std::endl;
      
      msg_Error() << std::setprecision(16)
                << "@@@ SHVIRT hasfsr=" << (HasFSR() ? 1 : 0)
                << " s=" << m_Sc.real() << " t=" << m_Tc.real()
                << " u=" << m_Uc.real() << " sQ=" << m_sQ
                << " deltI=(" << m_vertexI.real() << "," << m_vertexI.imag() << ")"
                << " deltF=(" << m_vertexF.real() << "," << m_vertexF.imag() << ")"
                << " BoxGGtu=(" << m_BoxGGtu.real() << "," << m_BoxGGtu.imag() << ")"
                << " BoxGZtu=(" << m_BoxGZtu.real() << "," << m_BoxGZtu.imag() << ")"
                << " BoxGGut=(" << m_BoxGGut.real() << "," << m_BoxGGut.imag() << ")"
                << " BoxGZut=(" << m_BoxGZut.real() << "," << m_BoxGZut.imag() << ")"
                << " rawGG=(" << m_rawboxGG.real() << "," << m_rawboxGG.imag() << ")"
                << " rawGZ=(" << m_rawboxGZ.real() << "," << m_rawboxGZ.imag() << ")"
                << " rawSub=(" << m_rawboxSub.real() << "," << m_rawboxSub.imag() << ")"
                << " coef=(" << m_rawcoef.real() << "," << m_rawcoef.imag() << ")"
                << " alpha=" << m_alpha << " MZ=" << m_MZ << " GZ=" << m_gZ
                << std::setprecision(10) << std::endl;
      msg_Error() << "@@@ SHPROP sp=" << m_sp
                << " propG=(" << m_propG.real() << "," << m_propG.imag() << ")"
                << " propZ=(" << m_propZ.real() << "," << m_propZ.imag() << ")"
                << std::endl;
      for (int j = 0; j <= 1; ++j)
        msg_Error() << "@@@ SHFF j=" << j
                  << " TC=(" << m_TC[j].real() << "," << m_TC[j].imag() << ")"
                  << " UC=(" << m_UC[j].real() << "," << m_UC[j].imag() << ")"
                  << std::endl;
      
      msg_Error() << "@@@ SHPHEL n=" << m_PhoHel.size()
                << " hel0=" << (m_PhoHel.empty()? 0 : m_PhoHel[0])
                << " nisr=" << m_isrphotons.size()
                << " nfsr=" << m_fsrphotons.size()
                << " nparts=" << m_nparts
                << std::endl;
      for (int h = 0; h <= 1; ++h) {
        const Complex sf(Sfactor(m_pceex[0], m_pceex[1], m_allphotons[0],
                                 h == 0 ? 1 : -1));
        msg_Error() << "@@@ SHSINI h=" << h
                  << " S=(" << sf.real() << "," << sf.imag() << ")"
                  << " |S|=" << std::abs(sf) << std::endl;
      }
      for (int j1 = 0; j1 <= 1; ++j1)
       for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
         for (int j4 = 0; j4 <= 1; ++j4) {
           const Complex tt(m_Tamp[j1][j2][j3][j4]);
           const Complex uu(m_Uamp[j1][j2][j3][j4]);
           const Complex bc(bref.m_A[Idx(j1,j2,j3,j4)]);
           if (std::abs(tt)==0. && std::abs(uu)==0. && std::abs(bc)==0.) continue;
           msg_Error() << "@@@ SHSPIN " << j1 << j2 << j3 << j4
                     << " TT=(" << tt.real() << "," << tt.imag() << ")"
                     << " UU=(" << uu.real() << "," << uu.imag() << ")"
                     << " Born=(" << bc.real() << "," << bc.imag() << ")"
                     << std::endl;
         }
    }

    
    msg_Error() << std::setprecision(10);
    for (int j1 = 0; j1 <= 1; ++j1)
      for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
          for (int j4 = 0; j4 <= 1; ++j4) {
            const Complex a0(m_AmpExpo0.m_A[Idx(j1,j2,j3,j4)]);
            const Complex a1(m_AmpExpo1.m_A[Idx(j1,j2,j3,j4)]);
            if (std::abs(a0) == 0. && std::abs(a1) == 0.) continue;
            const Complex bo(m_snapBorn.m_A[Idx(j1,j2,j3,j4)]);
            const Complex vi(m_snapVirt.m_A[Idx(j1,j2,j3,j4)]);
            const Complex re(m_snapReal.m_A[Idx(j1,j2,j3,j4)]);
            msg_Error() << "@@@ SHPART " << j1 << j2 << j3 << j4
                      << " born=(" << bo.real() << "," << bo.imag() << ")"
                      << " virt=(" << vi.real() << "," << vi.imag() << ")"
                      << " real=(" << re.real() << "," << re.imag() << ")"
                      << " virt/born=" << (std::abs(bo)>0.? std::abs(vi)/std::abs(bo):0.)
                      << " real/born=" << (std::abs(bo)>0.? std::abs(re)/std::abs(bo):0.)
                      << std::endl;
            msg_Error() << "@@@ SHAMP " << j1 << j2 << j3 << j4
                      << " A0=(" << a0.real() << "," << a0.imag() << ")"
                      << " A1=(" << a1.real() << "," << a1.imag() << ")"
                      << " A1/A0=(" << (std::abs(a0)>0.? (a1/a0).real() : 0.)
                      << "," << (std::abs(a0)>0.? (a1/a0).imag() : 0.) << ")"
                      << std::endl;
          }
    msg_Out() << "=== Sherpa CEEX point written to ceex_point.dat ===\n"
              << "  rho0 (O(alpha^0)) = " << m_result0 << "\n"
              << "  rho1 (O(alpha^1)) = " << m_result1 << "\n"
              << "  rho1/rho0 - 1     = " << (m_result1/m_result0 - 1.) << "\n"
              << "  compare against the KKMC harness's\n"
              << "    (RhoExp1 - RhoExp0)/RhoExp0 on the same point.\n"
              << "  sqrt(s) of this point = "
              << (m_momenta[0]+m_momenta[1]).Mass()
              << "   (momentum balance " << worst << ")\n"
              << "  m_sp (propagator s used)  = " << m_sp << "\n"
              << "  (P-k)^2 = s-channel momentum for ISR = "
              << (m_momenta[4]+m_momenta[5]).Abs2()
              << "   <-- CEEX eq.(one-photon): B(X) with X = P-k_1\n"
              << "  photon x = 2E/sqrt(s) = " << xgam
              << ",  E_gamma = " << m_allphotons[0][0] << " GeV\n";
    if (m_order == 2) {
      /*
        The tree-level O(alpha^2) partner of the KKMC harness's mode 5
        (KKTREE2 line: rho1/rho0 and rho2/rho0 with no virtual anywhere). The
        Sherpa side must be run without a CEEX-internal virtual
        (YFS: CEEX_Virtual: external, which composes outside m_result).
      */
      std::cerr<<std::setprecision(12)<<"@@@ SHTREE2 nphot="<<m_allphotons.size()
               <<" rho1/rho0="<<m_result1/m_result0
               <<" rho2/rho0="<<m_result2/m_result0
               <<" rho2/rho1="<<(m_result1 != 0. ? m_result2/m_result1 : 0.)
               <<" ownvirt="<<(CeexOwnVirtual() ? 1 : 0)
               <<" nisr="<<m_isrphotons.size();
      for (size_t i(0); i < m_allphotons.size(); ++i)
        std::cerr<<" x"<<i<<"="<<2.*m_allphotons[i][0]/(m_momenta[0]+m_momenta[1]).Mass();
      std::cerr<<std::endl;
    }
  }

  static const bool cxchk(ATOOLS::Settings::GetMainSettings()["CEEX"]["COMIX_CHECK"].Get<int>()!=0);
  if (cxchk) {
    ComixCalib c;
    if (CalibrateComixMap(c))
      std::cerr<<std::setprecision(10)
               <<"@@@ CEEXCMP sum_hand="<<c.sh<<" sum_comix="<<c.sc
               <<" comix_me2="<<c.me2
               <<" sc_over_me2="<<(c.me2!=0.? c.sc/c.me2 : 0.)
               <<" N="<<c.N
               <<" mask="<<c.mask<<" met="<<c.met<<" next="<<c.next
               <<" nlive="<<c.nlive
               <<" absratio=["<<c.rmin<<","<<c.rmax<<"]"
               <<" m_born="<<m_born
               <<" cth="<<(m_bornmomenta.size()>2?cos(m_bornmomenta[2].Theta()):0.)
               <<std::endl;
  }

  if (m_checkxs && m_bornmomenta.size() >= 4) {
    Vec4D_Vector bp(m_bornmomenta);
    Poincare bcms(bp[0] + bp[1]);
    for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);
    const double sp_save(m_sp);
    m_sp = (bp[2] + bp[3]).Abs2();
    MakeProp();
    MakePropT(bp);
    Amplitude b;
    BornAmplitude(bp, b);
    double sum(0.);
    for (int j1 = 0; j1 <= 1; ++j1)
      for (int j2 = 0; j2 <= 1; ++j2)
        for (int j3 = 0; j3 <= 1; ++j3)
          for (int j4 = 0; j4 <= 1; ++j4) {
            const Complex a(m_e * m_e * b.m_A[Idx(j1,j2,j3,j4)]);
            sum += std::real(a * conj(a));
          }
    m_sp = sp_save;
    MakeProp();
    MakePropT(m_pceex);
    AccumulateBornShape(cos(bp[2].Theta()), sum / 4.);
  }
}
