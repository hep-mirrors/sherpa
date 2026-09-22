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
  static const char *pin(getenv("CEEX_PIN_PHOTON_HEL"));
  static const int pinned(pin ? atoi(pin) : 0);
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
  m_svarQ = (m_pceex[2] + m_pceex[3]).Abs2();
  m_sQ    = m_svarQ;
  MakePropT(m_pceex);

  
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
    m_cfac = sProd * Complex(svarX/m_svarQ, 0.);
    MakeProp();
    // Electroweak form factors at THIS partition's scale, and the scattering
    // angle the WW/ZZ boxes depend on. No-op unless CEEX: WEAK is set.
    MakeEWFF(m_sp, cos(m_pceex[2].Theta()));
    MakeBoxMandelstams(PX);

    InfraredSubtractedME_0_0();
    InfraredSubtractedME_0_1();
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
        const double CKine((m_pceex[2] + m_pceex[3]
                            + m_allphotons[j]).Abs2() / m_svarQ);
        InfraredSubtractedME_1_0_FSR(m_allphotons[j], m_PhoHel[j],
                                     sProd, Sactu[j], CKine, 3 + 4*(int)j, (int)j);
      }
    }

    if (last == 1) break;
    PartitionPlus(last);
    if (last == 2) break;
  }
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

  static const double dumpxmin(getenv("SHERPA_CEEX_DUMP_XMIN") ?
                               atof(getenv("SHERPA_CEEX_DUMP_XMIN")) : 0.);
  static const size_t dumpn(getenv("SHERPA_CEEX_DUMP_NPHOT") ?
                            atoi(getenv("SHERPA_CEEX_DUMP_NPHOT")) : 1);
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

  static const bool cxchk(getenv("SHERPA_CEEX_COMIX")!=NULL);
  if (cxchk && m_bornmomenta.size() >= 4 && p_bornproc) {
    msg_Debugging()<<METHOD<<"(): entering Comix Born comparison\n";
    Vec4D_Vector bp(m_bornmomenta);
    Poincare bcms(bp[0] + bp[1]);
    for (size_t i(0); i < bp.size(); ++i) bcms.Boost(bp[i]);
    Amplitude cx, hand;
    const double sp_save(m_sp);
    m_sp = (bp[2] + bp[3]).Abs2();
    MakeProp();
    BornAmplitude(bp, hand);
    const bool ok(ComixBornAmplitude(bp, cx));
    m_sp = sp_save;
    MakeProp();
    if (ok) {
      double sh(0.), sc(0.), worst(0.);
      for (int a = 0; a <= 1; ++a)
        for (int b = 0; b <= 1; ++b)
          for (int c = 0; c <= 1; ++c)
            for (int d = 0; d <= 1; ++d) {
              const Complex H(m_e*m_e*hand.m_A[Idx(a,b,c,d)]), C(cx.m_A[Idx(a,b,c,d)]);
              sh += std::norm(H); sc += std::norm(C);
              const double den(std::abs(H) + std::abs(C));
              if (den > 0.) worst = Max(worst, std::abs(H - C)/den);
            }
      std::cerr<<"@@@ CEEXCMP sum_hand="<<sh<<" sum_comix="<<sc
               <<" ratio="<<(sc!=0.? sh/sc : 0.)
               <<" worst_elem_reldiff="<<worst
               <<" t="<<m_tinv<<" cth="<<cos(bp[2].Theta())<<std::endl;
      static int nrow(0);
      if (nrow++ < 3)
        for (int a = 0; a <= 1; ++a)
          for (int b = 0; b <= 1; ++b)
            for (int c = 0; c <= 1; ++c)
              for (int d = 0; d <= 1; ++d) {
                const Complex H(m_e*m_e*hand.m_A[Idx(a,b,c,d)]), C(cx.m_A[Idx(a,b,c,d)]);
                std::cerr<<"@@@ CEEXHEL "<<a<<b<<c<<d
                         <<" absH="<<std::abs(H)<<" absC="<<std::abs(C)
                         <<" H=("<<H.real()<<","<<H.imag()<<")"
                         <<" C=("<<C.real()<<","<<C.imag()<<")"<<std::endl;
              }
    }
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
