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

using namespace YFS;


bool Ceex_Base::BuildStages() {
  /*
    Stage 1 = initial, stage 0 = final: the index convention the ISR/FSR
    version used, kept so the odometer enumerates partitions in the same order.
    A 2 -> N core replaces the four pushes below with one stage per resonant
    propagator, and nothing downstream changes, because everything reads
    m_stagelegs.

    theta = -1 incoming, +1 outgoing, so w = Q*theta reproduces the weights
    Sfactor already applies implicitly: (+1,-1) for the incoming e-/e+ pair,
    and for an outgoing f/fbar pair the overall sign and magnitude that the
    explicit qratio = -q_f/q_e used to carry.
  */
  if (m_flavs.size() < 4) return false;
  // External legs first, in the order the stage legs will name them. Initial
  // from the Born momenta, final from the CEEX momenta. A 2 -> N core appends
  // its reconstructed resonances after these.
  m_stagemom.clear();
  for (size_t i(0); i < m_flavs.size(); ++i)
    m_stagemom.push_back(i < 2 ? m_bornmomenta[i] : m_pceex[i]);

  /*
    The FLAT decomposition: every initial leg in the production stage, every
    final leg in the decay stage. Valid for any N - the incoming pair is
    charge neutral and the final state is neutral by charge conservation - and
    at 2 -> 2 it is leg for leg the ISR/FSR split it replaces.

    This is deliberately NOT the resonance decomposition. Sub-staging a final
    state by its resonances (production, then one stage per resonance) is the
    G > 2 case and needs the resonance grouping; until then a 2 -> N process
    gets the flat scheme, which is a valid decomposition, just not the one
    that resums initial-final interference through a resonance.
  */
  m_stagelegs.assign(2, std::vector<StageLeg>());
  for (size_t i(0); i < m_flavs.size(); ++i) {
    const double w(m_flavs[i].Charge() * (i < 2 ? -1. : +1.));
    m_stagelegs[i < 2 ? 1 : 0].push_back({(int)i, w});
  }

  /*
    A stage factor is gauge invariant only if that stage is separately
    CHARGE-NEUTRAL: under eps -> eps + lambda k it varies by the partial sum
    sum_{i in g} Q_i theta_i, which the total charge conservation of the
    process does not make vanish. A decomposition failing this is rejected
    rather than silently producing a gauge-dependent answer.
  */
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
  m_stagereduces[1] = 1;   // only initial-stage photons reduce X
  return true;
}


void Ceex_Base::CalculateSfactors() {
  if (!BuildStages()) { m_Sfac.clear(); return; }
  m_Sfac.assign(m_nstages, std::vector<Complex>());
  for (size_t i(0); i < m_allphotons.size(); ++i)
    for (int g(0); g < m_nstages; ++g) {
      /*
        The stage's eikonal current, summed over its own legs against their
        weights w = Q*theta. Any number of legs; a charge-neutral pair
        reproduces the old single Sfactor() call up to the order the rounding
        happens, which is why the reference stream shifts at 1e-16 and not
        above it.

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

  /*
    Closure, eq. (6.4): the stage currents must sum to the TOTAL eikonal,

        sum_g s_g(k)  =  s(k)  =  sum_{all legs i} w_i * SfactorLeg(p_i,k)

    which is the statement that every emitter is counted exactly once. It is
    what catches a leg assigned to two stages, one left out, a wrong theta
    sign, or a resonance whose charge does not match its daughters -- an
    internal resonance leg has to appear with opposite weight in its production
    and its decay stage and cancel here.

    Charge neutrality per stage (4.4) is necessary but NOT sufficient for this:
    a decomposition can be neutral stage by stage and still double count.

    Gated on CHECK_XS so production pays nothing.
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


void Ceex_Base::PartitionStart(int &last) {
  const bool fsr(HasFSR());
  m_stage.assign(m_allphotons.size(), fsr ? 0 : 1);
  last = fsr ? 0 : 1;
}


void Ceex_Base::PartitionPlus(int &last) {
  const size_t n(m_stage.size());
  if (n == 0) { last = 2; return; }
  // Odometer in base m_nstages, least significant photon first. At
  // m_nstages = 2 this visits partitions in the same order as the original
  // ISR/FSR increment it replaces.
  size_t i(0);
  for (; i < n; ++i) {
    if (++m_stage[i] < m_nstages) break;
    m_stage[i] = 0;
  }
  if (i == n) { last = 2; return; }   // carried off the top: enumeration done
  // last = 1 marks the final partition, so the caller processes it and stops.
  bool atmax(true);
  for (size_t j(0); j < n; ++j)
    if (m_stage[j] != m_nstages - 1) { atmax = false; break; }
  last = atmax ? 1 : 0;
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
  m_beta10 = 0.0;
  m_beta01 = 0.0;
  m_beta00 = 0.0;
  BuildCeexMomenta();
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

  int last(0), nparts(0);
  
  static const size_t maxphot(Settings::GetMainSettings()["CEEX"]["MAX_PARTITION_PHOTONS"]
                              .SetDefault(100).Get<int>());
  if (m_allphotons.size() > maxphot) {
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): "<<m_allphotons.size()<<" photons would need "
                 <<"2^"<<m_allphotons.size()<<" partitions. Falling back to "
                 <<"the all-ISR partition for events above CEEX: "
                 <<"MAX_PARTITION_PHOTONS = "<<maxphot<<", which DROPS the "
                 <<"initial-final interference on them. Reported once."
                 <<std::endl;
    }
    m_stage.assign(m_allphotons.size(), 1);
    last = 1;
  } else {
    PartitionStart(last);
  }
  if (m_allphotons.size() > m_maxnphot) m_maxnphot = m_allphotons.size();
  for (;;) {
    ++nparts;
    Vec4D PX(m_pceex[0] + m_pceex[1]);
    Complex sProd(1., 0.);
    std::vector<Complex> Sactu(m_allphotons.size(), Complex(1., 0.));
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      Sactu[j] = m_Sfac[m_stage[j]][j];
      sProd   *= Sactu[j];
      if (m_stagereduces[m_stage[j]]) PX -= m_allphotons[j];
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
    m_cfac = sProd * Complex(svarY/m_svarQ, 0.);
    MakeProp();
    // Electroweak form factors at THIS partition's scale, and the scattering
    // angle the WW/ZZ boxes depend on. No-op unless CEEX: WEAK is set.
    MakeEWFF(m_sp, cos(m_pceex[m_if1].Theta()));
    MakeBoxMandelstams(PX);

    InfraredSubtractedME_0_0();
    // Only when a virtual was asked for and CEEX is its source; see
    // Ceex_Base::CeexOwnVirtual.
    if (CeexOwnVirtual()) InfraredSubtractedME_0_1();
    for (size_t j(0); j < m_allphotons.size(); ++j) {
      if (m_checkxs && m_stagereduces[m_stage[j]]) {
        CheckSoftFactor(m_allphotons[j]);
        CheckUVDiagonality(m_allphotons[j]);
      }
      if (m_stagereduces[m_stage[j]]) {
        InfraredSubtractedME_1_0(m_allphotons[j], m_PhoHel[j], sProd, Sactu[j],
                                 1 + 4*(int)j, (int)j);
      } else {
        // CKine = (q1+q2+k)^2/(q1+q2)^2, KKMC's svarX1/svarQ.
        const double CKine((m_pceex[m_if1] + m_pceex[m_if2]
                            + m_allphotons[j]).Abs2() / m_svarQ);
        InfraredSubtractedME_1_0_FSR(m_allphotons[j], m_PhoHel[j],
                                     sProd, Sactu[j], CKine, 3 + 4*(int)j, (int)j);
      }
    }

    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
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
 
  {
    const size_t n(m_allphotons.size());
    const int want(!HasFSR() || n > maxphot ? 1
                   : (n < 31 ? (1 << n) : nparts));
    if (nparts != want) {
      ++m_partbad;
      static bool warned(false);
      if (!warned) {
        warned = true;
        msg_Error()<<METHOD<<"(): partition counter enumerated "<<nparts
                   <<" partitions for "<<n<<" photons, expected "<<want
                   <<". Reported once."<<std::endl;
      }
    }
    ++m_partn;
  }

  // Comix-supplied O(alpha) real, if CEEX: COMIX_REAL asked for it. It
  // OVERWRITES m_AmpExpo1 for one-photon events, so it has to run after the
  // partition loop has finished accumulating and before MakeRho squares.
  ApplyComixReal();

  MakeRho();

  static const double dumpxmin(ATOOLS::Settings::GetMainSettings()["CEEX"]["DUMP_XMIN"].Get<double>());
  static const size_t dumpn(ATOOLS::Settings::GetMainSettings()["CEEX"]["DUMP_NPHOT"].Get<int>());
  const double xgam(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                    2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass() : 0.);
  if (m_checkxs && !m_ceexdumped && m_allphotons.size() == dumpn &&
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
              << "  rho1 (O(alpha^1)) = " << m_result << "\n"
              << "  rho1/rho0 - 1     = " << (m_result/m_result0 - 1.) << "\n"
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
