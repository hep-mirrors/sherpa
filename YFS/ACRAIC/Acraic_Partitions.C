#include "YFS/ACRAIC/Acraic_Base.H"
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


bool Acraic_Base::BuildStages() {
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
    m_stagemom.push_back(i < 2 ? m_bornmomenta[i] : m_pacraic[i]);

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
      for (size_t k(0); k < r.daughters.size(); ++k) pw += m_pacraic[r.daughters[k]];
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
    charge-neutral. A dipole group that is not (a same-charge pair from the
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
    The production stage: the beams, and every resonance as an outgoing
    leg. With no resonance this is the initial stage as before. The W terms
    cancel between here and the decay stages in the total eikonal, so the
    closure below holds with any common W momentum.
  */
  std::vector<StageLeg> ini;
  for (int i(0); i < 2; ++i) ini.push_back({i, -m_flavs[i].Charge()});
  for (size_t r(0); r < m_resonances.size(); ++r)
    ini.push_back({m_resonances[r].stagemomIndex, m_resonances[r].flav.Charge()});
  m_stagelegs.push_back(ini);
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
          // braced: msg_Debugging() is itself an if, and would take the else
          if (IsResonanceLeg(leg)) {
            msg_Debugging()<<ResonanceOf(leg).flav
                           <<(m_stagelegs[g][i].w*ResonanceOf(leg).flav.Charge() < 0.
                              ? "(in)" : "(out)");
          }
          else {
            msg_Debugging()<<leg;
          }
        }
        msg_Debugging()<<"}"<<(g==(size_t)m_initstage?"(production)":""); }
      msg_Debugging()<<"\n"; } }
  if (WStagesActive()) {
    static bool shown(false);
    if (!shown) { shown = true;
      msg_Info()<<"ACRAIC: W stages active for "<<m_flavs.size()-2
                <<"-fermion final state: "<<m_nstages<<" stages (production with "
                <<m_resonances.size()<<" W legs, one decay stage per W)."
                <<std::endl; }
  }
  return true;
}


bool Acraic_Base::IsResonanceLeg(int leg) const
{
  return leg >= (int)m_flavs.size() && leg < (int)m_stagemom.size();
}

const Acraic_Base::Resonance &Acraic_Base::ResonanceOf(int leg) const
{
  for (size_t r(0); r < m_resonances.size(); ++r)
    if (m_resonances[r].stagemomIndex == leg) return m_resonances[r];
  THROW(fatal_error, "stage leg "+ToString(leg)+" names no resonance");
}

std::vector<int> Acraic_Base::StageExternalLegs(int stage) const
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

size_t Acraic_Base::StageShiftMask(int stage) const
{
  size_t mask(0);
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return mask;
  if (stage == m_initstage) {
    for (size_t l(0); l < m_stagelegs[stage].size(); ++l) {
      const int leg(m_stagelegs[stage][l].leg);
      if (leg >= 0 && leg < 2) mask = AddLeg(mask, leg);
    }
    return mask;
  }
  const std::vector<int> legs(StageExternalLegs(stage));
  for (size_t l(0); l < legs.size(); ++l)
    if (legs[l] >= 0 && legs[l] < (int)m_flavs.size()) mask = AddLeg(mask, legs[l]);
  return mask;
}

void Acraic_Base::AddPlacementVetoes(const std::vector<int> &photons,
                                   PropShifts &sh) const
{
  if (m_stage.size() != m_allphotons.size()) return;
  // every placement of the explicit photons: which of them sit on stage g
  std::vector<std::vector<bool> > placements(1, std::vector<bool>(photons.size(), false));
  for (size_t q(0); q < photons.size(); ++q) {
    const size_t n(placements.size());
    for (size_t a(0); a < n; ++a) {
      std::vector<bool> on(placements[a]);
      on[q] = true;
      placements.push_back(on);
    }
  }
  const size_t first(m_flavs.size());   // Comix leg of the first explicit photon
  for (int g(0); g < m_nstages; ++g) {
    if (g == m_initstage) continue;
    const size_t stagelegs(StageShiftMask(g));
    if (stagelegs == 0) continue;
    std::vector<bool> actual(photons.size(), false);
    for (size_t q(0); q < photons.size(); ++q)
      actual[q] = (m_stage[photons[q]] == g);
    for (size_t a(0); a < placements.size(); ++a) {
      if (placements[a] == actual) continue;
      size_t legs(stagelegs);
      for (size_t q(0); q < photons.size(); ++q)
        if (placements[a][q]) legs = AddLeg(legs, first + q);
      sh.push_back(PHASIC::Process_Base::VetoEntry(legs));
    }
  }
}

std::string Acraic_Base::NoIFIRefusal()
{
  if (!HasFSR()) return "";   // nothing to interfere with
  if (m_nstages != 2)
    return "it needs exactly one initial and one final stage (e+e- -> f fbar), "
           "this process has "+ToString(m_nstages)+" stages";
  if (m_bhabha || m_tchannel != 0 || BornHasExchangeLine())
    return "the Born has a t-channel exchange line; the split by emitting "
           "stage is not gauge invariant there";
  if (WStagesActive()) return "W stages are active";
  if (ResonanceExponentiation())
    return "ACRAIC: RESONANCE_EXPONENTIATION exponentiates an initial-final "
           "logarithm";
  static const beta1legs::code legsmode(Settings::GetMainSettings()["ACRAIC"]
                                        ["BETA1_LEGS"].Get<beta1legs::code>());
  if (legsmode == beta1legs::balanced)
    return "ACRAIC: BETA1_LEGS balanced evaluates M_1 without propagator "
           "shifts, which carry the placement vetoes";
  if (m_hasvirtual && m_acraicvirtsrc != acraicvirt::acraic
      && m_acraicvirtsrc != acraicvirt::none)
    msg_Error()<<"YFS: IFI false with an external virtual: its initial-final "
               <<"boxes (IR-subtracted) stay in both columns. For no IFI at "
               <<"all use YFS: ACRAIC_Virtual: acraic (vertices only)."<<std::endl;
  msg_Info()<<"YFS: IFI false - no initial-final interference: incoherent "
            <<"partitions, M_1/M_2 split by emitting stage, no IF form factor, "
            <<"no IF boxes in ACRAIC's own virtual."<<std::endl;
  return "";
}

Vec4D Acraic_Base::StageSystemMomentum(int stage) const
{
  Vec4D q;
  if (stage < 0 || stage >= (int)m_stagelegs.size()) return q;
  const std::vector<int> legs(StageExternalLegs(stage));
  for (size_t l(0); l < legs.size(); ++l)
    if (legs[l] >= 0 && legs[l] < (int)m_pacraic.size()) q += m_pacraic[legs[l]];
  return q;
}

Vec4D Acraic_Base::StagePhotonSum(int stage, int iphot, int lphot) const
{
  Vec4D K;
  if (m_stage.size() != m_allphotons.size()) return K;
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if ((int)i != iphot && (int)i != lphot && m_stage[i] == stage)
      K += m_allphotons[i];
  return K;
}

Acraic_Base::LegMomentum
Acraic_Base::StageLegMomentum(const StageLeg &l, int iphot,
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
  if (m_weikonal == weikonal::partition && r.decayStage >= 0)
    m.pole += StagePhotonSum(r.decayStage, iphot);
  return m;
}

Complex Acraic_Base::StageCurrent(int stage, int iphot, const Vec4D &k, int hel,
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

void Acraic_Base::RecomputeResonanceSfactors()
{
  if (!WStagesActive() || m_weikonal == weikonal::daughters) return;
  if ((int)m_Sfac.size() != m_nstages) return;
  for (int g(0); g < m_nstages; ++g) {
    bool hasres(false);
    for (size_t l(0); l < m_stagelegs[g].size(); ++l)
      if (IsResonanceLeg(m_stagelegs[g][l].leg)) hasres = true;
    if (!hasres) continue;
    for (size_t i(0); i < m_allphotons.size() && i < m_Sfac[g].size(); ++i)
      m_Sfac[g][i] = StageCurrent(g, (int)i, m_allphotons[i], m_PhoHel[i], m_pacraic);
  }
}


void Acraic_Base::CalculateSfactors() {
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
*/
bool Acraic_Base::RepairMomentumBalance()
{
  static const bool on(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                      ["MOMENTUM_REPAIR"].Get<bool>());
  m_repairresid = -1.;
  if (!on || m_pacraic.size() < 3) return false;
  std::vector<Vec4D*> v;
  std::vector<double> m2;
  for (size_t i(2); i < m_pacraic.size(); ++i) {
    v.push_back(&m_pacraic[i]);
    m2.push_back(sqr(i < m_flavs.size() ? m_flavs[i].Mass() : 0.));
  }
  for (size_t i(0); i < m_isrphotons.size(); ++i) { v.push_back(&m_isrphotons[i]); m2.push_back(0.); }
  for (size_t i(0); i < m_fsrphotons.size(); ++i) { v.push_back(&m_fsrphotons[i]); m2.push_back(0.); }
  const size_t n(v.size());
  if (n < 2) return false;
  const Vec4D P(m_pacraic[0] + m_pacraic[1]);
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
    int piv[4] = {0, 1, 2, 3};
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
    }
    for (size_t i(0); i < n; ++i) {
      Vec4D &p(*v[i]);
      const double u[4] = {p[0], -p[1], -p[2], -p[3]};
      const double mu((0.5*r[i] - (u[0]*L[0] + u[1]*L[1] + u[2]*L[2] + u[3]*L[3]))*inv[i]);
      for (int a(0); a < 4; ++a) p[a] += L[a] + mu*u[a];
    }
  }
  return true;
}


void Acraic_Base::PartitionStart(int &last) {
  const bool fsr(HasFSR());
  static const double softcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                              ["SOFT_PARTITION_CUT"].Get<double>());
  // YFS: IFI false: no collapse - a collapsed photon carries the coherent
  // total eikonal, i.e. its soft initial-final interference
  const double cut(m_ifi ? softcut : 0.);
  const size_t n(m_allphotons.size());
  m_stage.assign(n, fsr ? 0 : m_initstage);
  m_fixedstage.assign(n, 0);
  size_t nfree(0);
  for (size_t i(0); i < n; ++i) {
    const double x(m_s > 0. ? 2.*m_allphotons[i][0]/sqrt(m_s) : 1.);
    if (fsr && x >= cut && (int)m_Sfac.size() == m_nstages) { ++nfree; continue; }
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


void Acraic_Base::PartitionPlus(int &last) {
  const size_t n(m_stage.size());
  if (n == 0) { last = 2; return; }
  // Odometer in base m_nstages over the enumerated photons only, least
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


void Acraic_Base::PrintKinematics(std::ostream &o) const
{
  o<<std::setprecision(6)<<"HEAVY legs:";
  for (size_t i(0); i < m_pacraic.size(); ++i)
    o<<" ["<<i<<"] "<<m_flavs[i]<<" E="<<m_pacraic[i][0]
     <<" th="<<m_pacraic[i].Theta()<<" m="<<m_pacraic[i].Mass();
  o<<std::endl;
  const double rs(m_s > 0. ? sqrt(m_s) : 1.);
  for (size_t j(0); j < m_allphotons.size(); ++j) {
    const Vec4D &k(m_allphotons[j]);
    double thmin(M_PI);
    for (size_t i(2); i < m_pacraic.size(); ++i)
      if (m_flavs[i].Charge() != 0.) {
        const double c(Vec3D(k)*Vec3D(m_pacraic[i])
                       /(Vec3D(k).Abs()*Vec3D(m_pacraic[i]).Abs()));
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


void Acraic_Base::MakePhotonHel() {
  static const int pinned(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["PIN_PHOTON_HEL"].Get<int>());
  m_PhoHel.clear();
  for (size_t i(0); i < m_allphotons.size(); ++i) {
    const int h(ran->Get() < 0.5 ? 1 : -1);
    m_PhoHel.push_back(pinned ? pinned : h);
  }
}


Acraic_Base::StageTableSwap::StageTableSwap(Acraic_Base &cb) :
  c(cb), legs(cb.m_stagelegs), mom(cb.m_stagemom), res(cb.m_resonances),
  sfac(cb.m_Sfac), stage(cb.m_stage), fixed(cb.m_fixedstage),
  reduces(cb.m_stagereduces), wgroups(cb.m_wgroups), groups(cb.m_groups),
  initstage(cb.m_initstage), nstages(cb.m_nstages) {}

Acraic_Base::StageTableSwap::~StageTableSwap()
{
  c.m_stagelegs = legs; c.m_stagemom = mom; c.m_resonances = res;
  c.m_Sfac = sfac; c.m_stage = stage; c.m_fixedstage = fixed;
  c.m_stagereduces = reduces; c.m_wgroups = wgroups; c.m_groups = groups;
  c.m_initstage = initstage; c.m_nstages = nstages;
}

/*
  Whether photon j is below the beta_1 threshold in EVERY partition. The cut
  (ComixInfraredSubtracted_1_0) is on y = 2k.S/S^2 with S the beams less the
  other initial-state photons, so the largest y is with all of them on the
  initial stage. A photon whose cut depends on the partition is left alone:
  the next-to-soft flux below is then neither added nor double counted.
*/
bool Acraic_Base::Beta1NeverEvaluated(size_t j) const
{
  static const double xcut(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                           ["BETA1_XCUT"].Get<double>());
  if (j >= m_allphotons.size() || m_pacraic.size() < 2 || WStagesActive())
    return false;
  Vec4D S(m_pacraic[0] + m_pacraic[1]);
  for (size_t i(0); i < m_allphotons.size(); ++i)
    if (i != j) S -= m_allphotons[i];
  const double S2(S.Abs2());
  if (!(S2 > 0.)) return false;
  return !(2.*(m_allphotons[j]*S)/S2 > xcut);
}

/*
  The next-to-soft flux of photon j on the stage this partition puts it on.
  beta_1 = M_1 - s beta_0 carries the photon's recoil on the Born; a photon
  without beta_1 has the eikonal alone, and the exact one-photon weight is
  the soft factor times (1 - x) for initial- and final-state photons alike.
  The initial stage: (X)^2/(X + k)^2 with X the beams less the partition's
  initial-state photons, i.e. s'/s. A final stage: (q + K)^2/(q + K - k)^2,
  this photon's share of the stage pseudoflux, which rho_1 does not carry
  (NO_PSEUDOFLUX rho0_only) while the crude does.
*/
double Acraic_Base::NextToSoftFlux(size_t j) const
{
  if (j >= m_allphotons.size() || j >= m_stage.size()) return 1.;
  const int g(m_stage[j]);
  const Vec4D &k(m_allphotons[j]);
  Vec4D with, without;
  if (g == m_initstage) {
    with = m_pacraic[0] + m_pacraic[1] - StagePhotonSum(g, -1);
    without = with + k;
  } else {
    if (g < 0 || g >= (int)m_stagelegs.size() || m_stagelegs[g].empty()) return 1.;
    without = StageSystemMomentum(g) + StagePhotonSum(g, (int)j);
    with = without + k;
  }
  const double a(with.Abs2()), b(without.Abs2());
  return (a > 0. && b > 0.) ? a/b : 1.;
}

double Acraic_Base::CrudeFromGenerator()
{
  static const bool crudeborn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                             ["CRUDE_BORN"].SetDefault(true).Get<bool>());
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  if (!crudeborn || !m_comixborn || !m_cxbalignok
      || m_pacraic.size() != m_flavs.size() || m_flavs.size() < 4) return -1.;
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
  const int nhel(Amplitude::NHel());
  double rho(0.);
  for (;;) {
    Vec4D PX(m_pacraic[0] + m_pacraic[1]);
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
      if (pfmode != pseudoflux::neither)
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
      const double crudered((pfmode == pseudoflux::neither ? 1. : pflux) * (m_s/X2)
                            * (crudefixed ? crudeprod : std::norm(sProd)));
      double rc(0.);
      for (int f(0); f < nhel; ++f)
        rc += crudered * std::norm(m_cxbalign.m_A[f]
                                   * Cred.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
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
  when more than ACRAIC: MAX_PARTITION_PHOTONS would be enumerated. The
  numerator partition sum and the generator-side crude (CrudeFromGenerator)
  both start from this, on their own stage tables.
*/
void Acraic_Base::PartitionStartCollapsed(int &last) {
  static const size_t maxphot(Settings::GetMainSettings()["ACRAIC"]["MAX_PARTITION_PHOTONS"]
                              .SetDefault(100).Get<int>());
  if (m_allphotons.size() > maxphot) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): "<<m_allphotons.size()<<" photons would need "
                 <<m_nstages<<"^"<<m_allphotons.size()<<" partitions. Above "
                 <<"ACRAIC: MAX_PARTITION_PHOTONS = "<<maxphot<<" the softest "
                 <<"photons are collapsed (fixed stage, total eikonal, as "
                 <<"below SOFT_PARTITION_CUT) until that many are left to "
                 <<"enumerate. Reported once."
                 <<(m_ifi ? "" : " With YFS: IFI false a collapsed photon "
                    "keeps its soft initial-final interference.")<<std::endl;
    }
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

void Acraic_Base::Calculate() {
  // Lab copies first: everything below boosts into the ACRAIC frame.
  m_allphotons_lab = m_isrphotons;
  m_allphotons_lab.insert(m_allphotons_lab.end(),
                          m_fsrphotons.begin(), m_fsrphotons.end());
  for (size_t i(0); i < m_isrphotons.size(); ++i) m_cms.Boost(m_isrphotons[i]);
  for (size_t i(0); i < m_fsrphotons.size(); ++i) m_cms.Boost(m_fsrphotons[i]);
  for (size_t i(0); i < m_bornmomenta.size(); ++i) m_cms.Boost(m_bornmomenta[i]);

  m_justdumped = false;
  m_rhocrud = 0.0;
  m_resexp.clear();
  m_resexpdv = 0.;
  m_resexpdh.clear();
  static const bool coh(Settings::GetMainSettings()["ACRAIC"]["IFI"]
                        .SetDefault(true).Get<bool>());
  m_ifi_coherent = coh || m_order == 2;
  /*
    YFS: IFI false: incoherent partitions at every order, each partition with
    its own placement of the photons in M_1 and M_2 (AddPlacementVetoes).
  */
  if (!m_ifi) {
    static const std::string why(NoIFIRefusal());
    if (!why.empty()) THROW(fatal_error, "YFS: IFI false is not available: "+why);
    m_ifi_coherent = false;
  }
  m_inc00 = m_inc11 = m_inc01 = 0.;
  m_inc22 = m_inc12 = m_inc02 = m_incbv = m_incbo = 0.;
  m_b1n = 0; m_b1min = m_b1max = m_b1sum = m_b1sq = 0.;  // beta_1 spread, per event
  m_beta10 = 0.0;
  m_beta01 = 0.0;
  m_beta00 = 0.0;
  BuildAcraicMomenta();
  RepairMomentumBalance();
  CheckBornNorm();
  ZerAmplit();
  m_allphotons = m_isrphotons;
  m_allphotons.insert(m_allphotons.end(),
                      m_fsrphotons.begin(), m_fsrphotons.end());
  MakePhotonHel();
  ProbeFrames();
  CalculateSfactors();
  /*
    CalculateSfactors clears m_Sfac when BuildStages rejects the decomposition
    (a stage carrying net charge, so gauge dependent). Bail here rather than
    index an empty table below: the event then contributes no ACRAIC weight,
    which YFS_Handler already treats as "ACRAIC produced nothing".
  */
  if ((int)m_Sfac.size() != m_nstages
      || (!m_allphotons.empty() && m_Sfac[0].size() != m_allphotons.size())) {
    m_nparts = 0;
    return;
  }
  BuildResonanceExponents();   // ACRAIC: RESONANCE_EXPONENTIATION
  m_spincache.resize(1 + 4*m_allphotons.size());
  m_spinvalid.assign(m_spincache.size(), 0);
  m_realphot.assign(m_allphotons.size(), Amplitude());
  m_realphotM1.assign(RealVirtualMode() == acraicrv::averaged ? m_allphotons.size() : 0,
                      Amplitude());
  SetBeta1Trace();
  /*
    The invariant of the outgoing FERMION pair, not of legs 2 and 3. With
    anything else in the final state those are not the same: for
    e+e- -> H l+ l- leg 2 is the Higgs, so this was (p_H + p_l-)^2 rather than
    the dilepton mass. m_svarQ is the denominator of the pseudo-flux factor
    m_cfac = sProd*(svarX/m_svarQ) and feeds the FSR kinematic term, so it is
    load bearing, and being wrong there is silent - it is a perfectly good
    invariant, just not the one the decay propagator sits at.
  */
  m_svarQ = (m_pacraic[m_if1] + m_pacraic[m_if2]).Abs2();
  m_sQ    = m_svarQ;
  MakePropT(m_pacraic);

  
  // One Comix/hand ratio per photon, before the partition loop: inside it the
  // soft weights are in scope but the photon's TOTAL amplitude is not.
  RunBeta1Scan();
  BuildComixBornAlignment();
  BuildResonanceBornDelta();   // after the alignment: it uses m_comixflip
  BuildComixPhotonRatios();
  // ACRAIC: TCHANNEL_REDUCED_BORN (tristate, see Acraic_Base::RegisterDefaults)
  ResolveReducedBorn();

  int last(0), nparts(0);
  PartitionStartCollapsed(last);
  if (m_allphotons.size() > m_maxnphot) m_maxnphot = m_allphotons.size();
  // ACRAIC: ORDER 2 - which photons can form beta_2 pairs (after the collapse)
  const bool dobeta2(m_order == 2 && PrepareBeta2());
  if (dobeta2) ++m_b2events;
  Amplitude snap0, snap1, snap2, snapbv, snapbo;
  // ACRAIC: NEXT_TO_SOFT_FLUX - which photons carry their next-to-soft flux
  static const bool nexttosoft(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                               ["NEXT_TO_SOFT_FLUX"].Get<bool>());
  std::vector<char> nextsoft(m_allphotons.size(), 0);
  if (nexttosoft)
    for (size_t j(0); j < m_allphotons.size(); ++j) nextsoft[j] = Beta1NeverEvaluated(j);
  for (;;) {
    ++nparts;
    if (!m_ifi_coherent) {
      snap0 = m_AmpExpo0; snap1 = m_AmpExpo1;
      snap2 = m_AmpBeta2; snapbv = m_AmpBornVirt; snapbo = m_AmpBornOnly;
    }
    Vec4D PX(m_pacraic[0] + m_pacraic[1]);
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
      // after the crude: it is the generator's density, not the amplitude
      if (nextsoft[j]) Sactu[j] *= sqrt(NextToSoftFlux(j));
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
    if (WStagesActive() && m_weikonal == weikonal::partition) {
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
        if (nextsoft[j]) Sactu[j] *= sqrt(NextToSoftFlux(j));
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

    Vec4D PY(PX);
    for (size_t i(2); i < m_pacraic.size(); ++i)
      if (i != m_if1 && i != m_if2) PY -= m_pacraic[i];
    const double svarY(PY.Abs2());
    m_svarY = svarY;          // the decay line's own partition-shifted scale
    m_PXvec = PX;             // and the momentum that realises m_sp
    
    static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                            ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
   
    m_pflux = 1.;
    if (pfmode != pseudoflux::neither)
      for (int g(0); g < m_nstages; ++g) {
        if (g == m_initstage || m_stagelegs[g].empty()) continue;
        // for a W decay stage q is the W's daughters: (P_W + K)^2/P_W^2
        const Vec4D q(StageSystemMomentum(g)), K(StagePhotonSum(g, -1));
        const double q2(q.Abs2());
        if (q2 > 0.) m_pflux *= (q + K).Abs2()/q2;
      }
    m_cfac = sProd * Complex(pfmode == pseudoflux::rho0_and_rho1 ? m_pflux : 1., 0.);   // rho_1 side
    MakeProp();
    ResonanceProbe(nparts);
    // Electroweak form factors at THIS partition's scale, and the scattering
    // angle the WW/ZZ boxes depend on. No-op unless ACRAIC: WEAK is set.
    MakeEWFF(m_sp, cos(m_pacraic[m_if1].Theta()));
    MakeBoxMandelstams(PX);

    AccumulateBeta1Partition();
    if (dobeta2 || m_treeprobe) ResetBeta1Record();
    InfraredSubtractedME_0_0();
    m_b1terms.clear();        // ACRAIC: TCHANNEL_MULTIPHOTON 2: this partition's
    // Only when a virtual was asked for and ACRAIC is its source; see
    // Acraic_Base::AcraicOwnVirtual.
    if (AcraicOwnVirtual()) InfraredSubtractedME_0_1();
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      if (m_checkxs && m_stagereduces[m_stage[j]]) {
        CheckSoftFactor(m_allphotons[j]);
        CheckUVDiagonality(m_allphotons[j]);
      }
      /*
        beta_1: Comix's one-photon amplitude with the eikonal subtracted off
        it.
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
    // ACRAIC: TCHANNEL_MULTIPHOTON 2: the factorised remainder of the beta_1
    // terms just added (empty, so a no-op, in every other mode)
    AddFactorisedRemainder();
    /*
      ACRAIC: ORDER 2 - the tree residuals beta_m, 2 <= m <= TREE_BETA_MAX, for
      every eligible photon set (Acraic_Beta2.C), after the beta_1 loop of the
      SAME partition, whose M_1 and eikonals they subtract.
    */
    if (dobeta2 || m_treeprobe) AddTreeBetas(dobeta2);
    EmulateKKMCFlux();

    if (!m_ifi_coherent)   // this partition's own contribution, squared
      for (int f(0); f < Amplitude::NHel(); ++f) {
        const Complex d0(m_AmpExpo0.m_A[f] - snap0.m_A[f]);
        const Complex d1(m_AmpExpo1.m_A[f] - snap1.m_A[f]);
        m_inc00 += std::norm(d0); m_inc11 += std::norm(d1);
        m_inc01 += std::real(std::conj(d0)*d1);
        // ORDER 2 (A_2 = A_1 + beta_2) and the Born+virtual, YFS: IFI false
        const Complex d2(d1 + m_AmpBeta2.m_A[f] - snap2.m_A[f]);
        const Complex dbv(m_AmpBornVirt.m_A[f] - snapbv.m_A[f]);
        m_inc22 += std::norm(d2);
        m_inc12 += std::real(std::conj(d1)*d2);
        m_inc02 += std::real(std::conj(d0)*d2);
        m_incbv += std::norm(dbv);
        m_incbo += std::norm(m_AmpBornOnly.m_A[f] - snapbo.m_A[f]);
      }
    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
  /*
    ACRAIC: CRUDE_FROM_GENERATOR - the crude on the generator's own stages.
    compare only prints it next to the per-partition one (the gate: the two
    must agree wherever ACRAIC's stages are the generator's dipoles).
  */
  const crudegen::code cg(m_crudegen != crudegen::automatic ? m_crudegen
                          : (WStagesActive() ? crudegen::on : crudegen::off));
  if (cg != crudegen::off) {
    const double rg(CrudeFromGenerator());
    if (cg == crudegen::compare) {
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
        msg_Error()<<METHOD<<"(): ACRAIC: CRUDE_FROM_GENERATOR asked for but the "
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
  CheckClosureOnePhoton();
  ReportBeta1Partitions();
  CheckBeta1Ratio();
  CheckVirtualPartitions(nparts);
  m_nparts = nparts;
  if (m_checkxs && nparts > 1 && m_pzmin > 0.) {
    const size_t n(Min(m_allphotons.size(), size_t(15)));
    m_pzspread[n] = Max(m_pzspread[n], m_pzmax/m_pzmin);
    m_pzsum[n] += m_pzmax/m_pzmin;
    ++m_pzn[n];
    const double ecm((m_momenta[0]+m_momenta[1]).Mass());
    for (int c(0); c < 4; ++c) {
      const double xcut(m_pzxcut[c]);
      Vec4D Xhard(m_pacraic[0] + m_pacraic[1]), Xall(Xhard);
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
      m_pzsoftsav[c] += double(NumLegSubsets<unsigned>(Min(nsoft, size_t(20))));
    }
  }
 
  CheckPartitionCount(nparts);

  // Comix-supplied O(alpha) real, if ACRAIC: COMIX_REAL asked for it. It
  // OVERWRITES m_AmpExpo1 for one-photon events, so it has to run after the
  // partition loop has finished accumulating and before MakeRho squares.
  ApplyComixReal();

  MakeRho();
  if (m_order == 2) {
    static const bool b2c(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["BETA2_CLOSURE"].Get<int>() != 0);
    static const bool b2s(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["BETA2_SOFT_TEST"].Get<int>() != 0);
    if (b2c) Beta2Closure();
    if (b2s) Beta2SoftLimit();
    if (m_treeprobe) TreeBetaClosure();   // ACRAIC: TREE_BETA_PROBE
  }
  RunBeta1Closure();
  RunSoftLimitTest();
  RunSoftNormCheck();

  static const double dumpxmin(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["DUMP_XMIN"].Get<double>());
  static const size_t dumpn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["DUMP_NPHOT"].Get<int>());
  const double xgam(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                    2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass() : 0.);
  // DUMP_XMIN_EACH: every photon above this x (a multi-hard-photon point)
  static const double dumpxeach(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["DUMP_XMIN_EACH"].Get<double>());
  bool dumpeach(true);
  if (dumpxeach > 0. && m_momenta.size() >= 2)
    for (size_t i(0); i < m_allphotons.size(); ++i)
      if (!(2.*m_allphotons[i][0]/(m_momenta[0]+m_momenta[1]).Mass() > dumpxeach))
        dumpeach = false;
  if (m_checkxs && !m_acraicdumped && m_allphotons.size() == dumpn && dumpeach &&
      xgam > dumpxmin && m_momenta.size() >= 6 && m_result0 != 0.) {
    
    Vec4D bal(m_momenta[0] + m_momenta[1] - m_momenta[4] - m_momenta[5]);
    for (size_t i(0); i < m_allphotons.size(); ++i) bal -= m_allphotons[i];
    const double worst(Max(Max(dabs(bal[0]), dabs(bal[1])),
                           Max(dabs(bal[2]), dabs(bal[3]))));
    if (worst > 1e-6) {
      msg_Error()<<METHOD<<"(): ACRAIC dump skipped, momentum imbalance "
                 <<worst<<std::endl;
      return;
    }
    m_acraicdumped = true;
    m_justdumped = true;
    std::ofstream pt("acraic_point.dat");
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
    
    Amplitude bref;
    BornAmplitude(m_pacraic, bref);
    msg_Error() << std::setprecision(10);
    msg_Error() << "@@@ SHMASS nflav=" << m_flavs.size()
              << " m1=" << (m_flavs.size()>0? m_flavs[0].Mass():-1.)
              << " m2=" << (m_flavs.size()>1? m_flavs[1].Mass():-1.)
              << " m3=" << (m_flavs.size()>2? m_flavs[2].Mass():-1.)
              << " m4=" << (m_flavs.size()>3? m_flavs[3].Mass():-1.)
              << " mass_I=" << m_mass_I << " mass_F=" << m_mass_F
              << std::endl;

    msg_Error() << std::setprecision(16) << "@@@ SHPMASS";
    for (size_t i = 0; i < m_pacraic.size(); ++i)
      msg_Error() << " p" << i << "=" << m_pacraic[i].Mass();
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
      const Complex sf(Sfactor(m_pacraic[0], m_pacraic[1], m_allphotons[0],
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
    msg_Out() << "=== Sherpa ACRAIC point written to acraic_point.dat ===\n"
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
              << "   <-- ACRAIC eq.(one-photon): B(X) with X = P-k_1\n"
              << "  photon x = 2E/sqrt(s) = " << xgam
              << ",  E_gamma = " << m_allphotons[0][0] << " GeV\n";
    if (m_order == 2) {
      /*
        The tree-level O(alpha^2) partner of the KKMC harness's mode 5
        (KKTREE2 line: rho1/rho0 and rho2/rho0 with no virtual anywhere). The
        Sherpa side must be run without a ACRAIC-internal virtual
        (YFS: ACRAIC_Virtual: external, which composes outside m_result).
      */
      std::cerr<<std::setprecision(12)<<"@@@ SHTREE2 nphot="<<m_allphotons.size()
               <<" rho1/rho0="<<m_result1/m_result0
               <<" rho2/rho0="<<m_result2/m_result0
               <<" rho2/rho1="<<(m_result1 != 0. ? m_result2/m_result1 : 0.)
               <<" ownvirt="<<(AcraicOwnVirtual() ? 1 : 0)
               <<" nisr="<<m_isrphotons.size();
      for (size_t i(0); i < m_allphotons.size(); ++i)
        std::cerr<<" x"<<i<<"="<<2.*m_allphotons[i][0]/(m_momenta[0]+m_momenta[1]).Mass();
      std::cerr<<std::endl;
    }
  }

  static const bool cxchk(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["COMIX_CHECK"].Get<int>()!=0);
  if (cxchk) {
    ComixCalib c;
    if (CalibrateComixMap(c))
      std::cerr<<std::setprecision(10)
               <<"@@@ ACRAICCMP sum_hand="<<c.sh<<" sum_comix="<<c.sc
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
    MakePropT(m_pacraic);
    AccumulateBornShape(cos(bp[2].Theta()), sum / 4.);
  }
}

void Acraic_Base::ResolveReducedBorn()
{
  static const tristate::code rb(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                 ["TCHANNEL_REDUCED_BORN"].Get<tristate::code>());
  int nrad(0);
  for (size_t g(0); g < m_stagelegs.size(); ++g)
    if (!m_stagelegs[g].empty()) ++nrad;
  m_redborn = rb == tristate::on || (rb == tristate::automatic && m_comixborn && nrad == 1
                         && BornHasExchangeLine());
  // the reduced point is LegsAt's 2 -> 2 construction; with W stages the
  // partition Born stays the shifted physical-spinor one
  if (WStagesActive() && m_redborn) {
    static bool warned(false);
    if (!warned) { warned = true;
      msg_Error()<<METHOD<<"(): ACRAIC: TCHANNEL_REDUCED_BORN is not defined "
                 <<"with W stages; ignored."<<std::endl; }
    m_redborn = false;
  }
}
