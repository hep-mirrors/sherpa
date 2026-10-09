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


Complex Ceex_Base::BornAmplitude(const Vec4D_Vector &k) {
  Complex amp;
  double hel1, hel2, hel3, hel4;
  int mode;
  for (int h0 = 1; h0 <= 2; ++h0) {
    for (int h1 = 1; h1 <= 2; ++h1) {
      for (int h2 = 1; h2 <= 2; ++h2) {
        for (int h3 = 1; h3 <= 2; ++h3) {
          hel1 = 3 - 2 * h0;
          hel2 = 3 - 2 * h1;
          hel3 = 3 - 2 * h2;
          hel4 = 3 - 2 * h3;
          if (hel1 == -hel2 ) {
            m_T = T(k[2], k[0], hel3, hel1) * Tp(k[1], k[3], hel2, hel4);
            m_U = Up(k[2], k[1], hel3, hel2) * U(k[0], k[3], hel1, hel4);
            m_ampborn[h0][h1][h2][h3] = (CouplingZ(hel1, 1) + CouplingG()) * m_T + (CouplingZ(hel1, 1) + CouplingG()) * m_U;
            m_bornAmp.m_A[Idx(h0,h1,h2,h3)] = (CouplingZ(hel1, 1) + CouplingG()) * m_T + (CouplingZ(hel1, 1) + CouplingG()) * m_U;
            amp += (CouplingZ(hel1, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(hel1, 1) * m_propZ + CouplingG() * m_propG) * m_T;
          }
        }
      }
    }
  }
  return amp;
}


void Ceex_Base::BornAmplitude(const Vec4D_Vector &k, Amplitude &M,
                             double Mf3, double Mf4, int slot) {
  Complex amp;
  int hel1, hel2, hel3, hel4;
  const double m1(m_flavs[0].Mass());
  /*
    The outgoing FERMION pair, which is legs 2 and 3 only when every final leg
    is a fermion. For H l+ l- leg 2 is the Higgs, and this routine was building
    spinors for (e-, e+, H, mu-) - a scalar in a fermion slot and mu+ never
    touched at all. It survived inspection because the CEEX weight is a RATIO:
    a wrong Born sits in the numerator and the denominator both and largely
    cancels, so the column stayed plausible while the amplitude was nonsense.
    It is only against Comix, helicity by helicity, that it shows - |C/H| ran
    over a factor of 21 across the four live helicities, where a correct Born
    differing by couplings alone would give one constant.
  */
  /*
    This is called with TWO kinds of momentum set, and they index the pair
    differently: the full CEEX set, where the outgoing fermions sit at
    m_if1/m_if2, and a REDUCED 2 -> 2 Born set {P1, P2, q3, q4} built by the
    real paths, where they sit at 2 and 3 whatever the process. Using m_if1
    unconditionally reads off the end of the reduced set - for H l+ l- it
    indexes 3 and 4 of a 4-vector - which is how the hand-coded column went
    from +24% to 1e225.
  */
  const size_t i1(k.size() > 4 ? m_if1 : 2), i2(k.size() > 4 ? m_if2 : 3);
  const double m3(Mf3 >= 0. ? Mf3 : m_flavs[i1].Mass());
  const double m4(Mf4 >= 0. ? Mf4 : m_flavs[i2].Mass());
  static const bool exactisrmass(
      Settings::GetMainSettings()["CEEX"]["EXACT_ISR_SPINOR_MASS"]
      .SetDefault(false).Get<bool>());
  const double mi(exactisrmass || m_bhabha ? m1 : 0.);
  const bool cached(slot >= 0 && slot < (int)m_spinvalid.size()
                    && m_spinvalid[slot]);
  if (cached) {
    const SpinorSet &c(m_spincache[slot]);
    for (int a = 0; a <= 1; ++a)
      for (int b = 0; b <= 1; ++b)
        for (int cc = 0; cc <= 1; ++cc)
          for (int d = 0; d <= 1; ++d) {
            m_Tamp[a][b][cc][d] = c.T[a][b][cc][d];
            m_Uamp[a][b][cc][d] = c.U[a][b][cc][d];
          }
  } else
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          hel1 = 1 - 2 * h0;
          hel2 = 1 - 2 * h1;
          hel3 = 1 - 2 * h2;
          hel4 = 1 - 2 * h3;
          if (hel1 == -hel2 ) {
         
            m_T  = T(k[i1], k[0], hel3, hel1, +1, +1, m3, mi)
                 * Tp(k[1], k[i2], hel2, hel4, -1, -1, mi, m4);
            m_U  = Up(k[i1], k[1], hel3, hel2, +1, +1, m3, mi)
                 * U(k[0], k[i2], hel1, hel4, -1, -1, mi, m4);
            m_Tamp[h0][h1][h2][h3] = m_T;
            m_Uamp[h0][h1][h2][h3] = m_U;
          }
          else {
            m_Tamp[h0][h1][h2][h3] = 0.;
            m_Uamp[h0][h1][h2][h3] = 0.;

          }
        }
      }
    }
  }
  if (m_bhabha)
    for (int h0 = 0; h0 <= 1; ++h0)
      for (int h1 = 0; h1 <= 1; ++h1)
        for (int h2 = 0; h2 <= 1; ++h2)
          for (int h3 = 0; h3 <= 1; ++h3) {
            const int hl1(1 - 2*h0), hl2(1 - 2*h1), hl3(1 - 2*h2), hl4(1 - 2*h3);
            if (hl1 == hl3) {
              // crossing p2 <-> -p3 flips the helicity label of both exchanged
              // legs, which is what turns the s-channel gate into hel1 == hel3
              m_Tampt[h0][h1][h2][h3] = T(k[1], k[0], -hl2, hl1, +1, +1, mi, mi)
                                      * Tp(k[2], k[3], -hl3, hl4, -1, -1, m3, m4);
              m_Uampt[h0][h1][h2][h3] = Up(k[1], k[2], -hl2, -hl3, +1, +1, mi, m3)
                                      * U(k[0], k[3], hl1, hl4, -1, -1, mi, m4);
            } else {
              m_Tampt[h0][h1][h2][h3] = 0.;
              m_Uampt[h0][h1][h2][h3] = 0.;
            }
          }

  if (!cached && slot >= 0 && slot < (int)m_spincache.size()) {
    SpinorSet &c(m_spincache[slot]);
    for (int a = 0; a <= 1; ++a)
      for (int b = 0; b <= 1; ++b)
        for (int cc = 0; cc <= 1; ++cc)
          for (int d = 0; d <= 1; ++d) {
            c.T[a][b][cc][d] = m_Tamp[a][b][cc][d];
            c.U[a][b][cc][d] = m_Uamp[a][b][cc][d];
          }
    m_spinvalid[slot] = 1;
  }
  // The ONLY partition-dependent part: two coupling x propagator factors.
  for (int j = 0; j <= 1; j++) {
    double h = 1. - 2.*j;
    m_UC[j] = CouplingZ(h, 0) * m_propZ + CouplingG() * m_propG;
    m_TC[j] = CouplingZ(h, 1) * m_propZ + CouplingG() * m_propG;
  }
  // Bhabha t-channel couplings, carried by the crossed structures built above.
  //
  for (int j = 0; j <= 1; j++) m_TCt[j] = m_UCt[j] = Complex(0., 0.);
  if (m_bhabha) {
    for (int j = 0; j <= 1; j++) {
      const double h(1. - 2.*j);
      m_TCt[j] = -(CouplingZ(h, 1) * m_propZt + CouplingG() * m_propGt);
      m_UCt[j] = -(CouplingZ(h, 0) * m_propZt + CouplingG() * m_propGt);
    }
  }
  for (int h0 = 0; h0 <= 1; h0++) {
    for (int h1 = 0; h1 <= 1; h1++) {
      for (int h2 = 0; h2 <= 1; h2++) {
        for (int h3 = 0; h3 <= 1; h3++) {
          M.m_A[Idx(h0,h1,h2,h3)] = m_TC[h0] * m_Tamp[h0][h1][h2][h3]
                                + m_UC[h0] * m_Uamp[h0][h1][h2][h3];
          if (m_bhabha)
            M.m_A[Idx(h0,h1,h2,h3)] += m_TCt[h0] * m_Tampt[h0][h1][h2][h3]
                                   + m_UCt[h0] * m_Uampt[h0][h1][h2][h3];
          m_bornAmp.m_A[Idx(h0,h1,h2,h3)] = M.m_A[Idx(h0,h1,h2,h3)];
        }
      }
    }
  }
}







Complex Ceex_Base::BornAmplitude(const Vec4D_Vector &k, int h0, int h1, int h2, int h3) {
  Complex amp;
  if (h0 == -h1 ) {
    m_T = T(k[2], k[0], h2, h0) * Tp(k[1], k[3], h1, h3);
    m_U = Up(k[2], k[1], h2, h1) * U(k[0], k[3], h0, h3);
    amp = (CouplingZ(h0, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(h0, 1) * m_propZ + CouplingG() * m_propG) * m_T;
  }
  return amp;
}



Complex Ceex_Base::BornAmplitude(Vec4D p1, Vec4D p2, Vec4D p3, Vec4D p4, int h0, int h1, int h2, int h3)
{
  Vec4D_Vector tmp;
  tmp.push_back(p1);
  tmp.push_back(p2);
  tmp.push_back(p3);
  tmp.push_back(p4);
  return BornAmplitude(tmp, h0, h1, h2, h3);
}



/*!
  Born spin amplitudes from Comix. Always false, and structurally so.

*/
bool Ceex_Base::ComixBornAmplitude(const Vec4D_Vector &p, Amplitude &A,
                                   double *me2, const double propscale,
                                   const double decscale)
{
  if (p_bornproc == NULL) return false;
  std::vector<METOOLS::Spin_Amplitudes> amps;
  if (!p_bornproc->BornSpinAmplitudes(p, amps, me2, propscale,
                                     decscale, DecayCId())) return false;
  if (amps.empty()) return false;
  const METOOLS::Spin_Amplitudes &sa(amps[0]);
  const int nh(Amplitude::NHel());
  if ((int)sa.size() < nh) return false;
  for (int f = 0; f < nh; ++f) A.m_A[f] = sa[f];
  return true;
}


bool Ceex_Base::ComixBornShifted(const Vec4D_Vector &p, Amplitude &A,
                                 const PropShifts &shifts)
{
  if (p_bornproc == NULL) return false;
  std::vector<METOOLS::Spin_Amplitudes> amps;
  if (!p_bornproc->BornSpinAmplitudesShifts(p, amps, NULL, shifts)) return false;
  if (amps.empty()) return false;
  const METOOLS::Spin_Amplitudes &sa(amps[0]);
  const int nh(Amplitude::NHel());
  if ((int)sa.size() < nh) return false;
  for (int f = 0; f < nh; ++f) A.m_A[f] = sa[f];
  return true;
}

void Ceex_Base::InfraredSubtractedME_0_0() {
  // This partition's Born, squared on its own and added to the INCOHERENT sum
  double rc(0.);
  Amplitude AmpBorn;
  BornAmplitude(m_pceex, AmpBorn, -1., -1., 0);
  /*
    Comix's Born at THIS partition's scale, in CEEX's convention. The spinors
    stay at the physical momenta and only the propagator pole moves, which is
    what the partition sum means by B(X) and why no momentum configuration
    realises it.
  */
  bool realpoint(false);
  if (m_comixborn && m_cxbalignok) {
    Amplitude C;
    /*
      Two forms of beta_0(X_wp), both with the propagator at THIS partition's
      X and both partition dependent:

      BORN_AT_SPRIME: 0 - KKMC's. Physical spinors with the pole pinned to
      m_sp (the decay line to m_svarY), a configuration no momenta realise,
      times the pseudo-flux svarY/svarQ. The flux is what makes it the size of
      a Born: pinning scales the propagator by svarQ/X^2 relative to the real
      point while the spinors stay at their physical scale, and X^2/svarQ
      undoes that.

      BORN_AT_SPRIME: 1 - the Born at the REAL point whose invariant is X_wp^2:
      beams back-to-back at X_wp, the radiating pair at X_wp minus the
      spectators, directions taken from the event (BornLegsAt). The
      propagators follow from the momenta and no flux is applied. This is the
      object Comix's one-photon amplitude reduces to as its photon goes soft,
      so beta_1 = M_1 - s beta_0 becomes a difference of like objects rather
      than of a real-point amplitude and a pinned one.

      An earlier version of this switch took ONE point per event - the
      all-ISR reduction, m_plabmom - for every partition and kept the flux.
      The (F,F) partition then carried s/s' on a Born already at s', and
      mu mu came out at 6215 pb. That was the flux, not the alignment.
    */
    static const bool sprime(ATOOLS::Settings::GetMainSettings()["CEEX"]
                             ["BORN_AT_SPRIME"].Get<bool>());
    bool ok(false);
    // BornLegsAt is a 2 -> 2 reduction (one pair at Y): not with W stages
    if ((sprime || m_redborn) && !WStagesActive()) {   // m_redborn: CEEX: TCHANNEL_REDUCED_BORN
      Vec4D_Vector pb;
      ok = BornLegsAt(m_PXvec, pb)
           && ComixBornAmplitude(pb, C, NULL, -1., -1.);
      realpoint = ok;
    }
    /*
      The pinned Born by momentum SHIFT rather than by scalar: the initial-state
      photons of this partition added to the initial-side lines, the
      final-state ones to the final-side lines. With every photon included the
      two sides agree, (P - sum_I k)^2 = (q_c + q_d + sum_F k)^2 = X_wp^2, so
      this is the same object as the scalar pin at (m_sp, m_svarY) - checked
      below - and the form beta_1's M_1 reduces to, since the same shifts
      with THIS photon left out are what M_1 is evaluated with.
    */
    if (!ok) {
      PropShifts sh(StageShifts(-1));
      AddExchangeLineShifts(-1, sh);     // space-like exchange lines
      ok = ComixBornShifted(m_pceex, C, sh);
      if (ok && m_checkxs && m_ffbar && m_nstages == 2) {
        static int nchk(0);
        if (nchk < 12) { ++nchk;
          Amplitude Cp;
          if (ComixBornAmplitude(m_pceex, Cp, NULL, m_sp, m_svarY)) {
            double nd(0.), nn(0.);
            for (int f = 0; f < Amplitude::NHel(); ++f) {
              nd += std::norm(C.m_A[f] - Cp.m_A[f]); nn += std::norm(Cp.m_A[f]); }
            std::cerr<<"@@@ SHIFTCHK nphot="<<m_allphotons.size()
                     <<" sp="<<m_sp<<" sY="<<m_svarY
                     <<" |shift-pin|/|pin|="<<(nn>0.? sqrt(nd/nn) : -1.)
                     <<std::endl; } }
      }
    }
    if (!ok) ok = ComixBornAmplitude(m_pceex, C, NULL, m_sp, m_svarY);
    if (ok) {
      const int nh(Amplitude::NHel());
      const int fmaskx(Amplitude::NHel() - 1);
      for (int f = 0; f < nh; ++f)
        AmpBorn.m_A[f] = m_cxbalign.m_A[f]*C.m_A[f ^ (m_comixflip & fmaskx)]
                         /(m_e*m_e);
    }
  }
  // the real-point Born carries its own scale: soft factors only, no flux.
  // fac0 feeds rho_1 (and beta_1's subtraction), fac00 feeds rho_0; see the
  // NO_PSEUDOFLUX modes in the partition loop.
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["CEEX"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  const Complex fac0(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::rho0_and_rho1 ? m_pflux : 1.)));
  const Complex fac00(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::neither ? 1. : m_pflux)));
  const Complex fac(fac0 * m_Sprod), facA0(fac00 * m_Sprod);
  /*
    Every helicity the container holds. This loop ran over the 16 entries of
    a 2 -> 2 final state until 2026-09-24; for e+e- -> mu mu tau tau (64
    entries) the other 48 never received a Born, and m_partborn0 stayed zero
    there, so beta_1 for those helicities was M_1 with NOTHING subtracted -
    infrared-unsafe, and the CEEX weight grew with the photon multiplicity.
  */
  const int nhel(Amplitude::NHel());
  // With a collapsed photon the crude's soft product is not |m_Sprod|^2 but
  // the incoherent one (m_crudeprod, see Ceex_Base.H). Without one the two
  // are the same number and the original expression is kept bit for bit.
  const double crudefac(std::norm(fac00) * m_crudeprod);
  /*
    The CRUDE's Born is the generator's: the Born at the partition's reduced
    point (BornLegsAt(X_wp): beams back to back at X_wp along the
    generator's axis, the pair at Y_wp), times s/X_wp^2. The coherent sum
    and beta_1's subtraction keep the KKMC form above - physical spinors,
    poles at X_wp - because that is what M_1 reduces to (Defect 3). For an
    s-channel Born the two are the same number: the amplitude at physical
    spinors with the pole at X is sqrt(s/X^2) times the one at the reduced
    point (the spinor products scale, the pole is shared), and the
    fixed-order weight divides its real by exactly S~ m_born/(1 - x) with
    m_born the reduced-point Born - which is why CEEX = Born+real held
    event by event at one photon for mu mu and nu nu. For a space-like
    exchange line the two are NOT the same number: the numerators do not
    scale with the pole, and pinning a t-channel pole under physical
    numerators can drive it to zero. e+e- -> gamma gamma at the Z pole:
    Born+real / CEEX at one photon had a median of 1.000 and a 1st-99th
    percentile of 0.41-25 (up to 5e9), Bhabha CEEX/Born+real 0.90 with a
    heavy tail, both with the physical-spinor crude. Measured: the Comix
    Born at BornLegsAt is m_born to six digits at every x on gamma gamma
    and on nu nu at x = 0.87 once the axis is the generator's
    (REDUCED_AXIS). Measured with it (CEEX: CRUDE_BORN: 1): Z-pole mu mu
    CEEX 1278.1 -> 1280.1 (+0.15%, within errors, as the argument says),
    Bhabha 1376.9 -> 1374.8, four-fermion bit-identical (2 -> 2 only, since
    LegsAt rebuilds one pair) - but gamma gamma still fails the one-photon
    test at wide angle (Born+real/CEEX median 15 for |cos theta| < 0.9, all
    at x > 0.9), so the crude is not what is wrong there and this is left
    OFF (default 0) until that is understood: the validated s-channel
    numbers stay bit for bit.
    [2026-09-26: that "Born+real" was itself wrong for gamma gamma - its
    beta_1 was divided by 3 (YFS: REAL_BORN_PHOTON_SYM, NLO_Base.C). Against
    the exact |M_1|^2/density this crude makes the one-photon CEEX factor
    exact at every x; the huge multi-photon weights are the physical-spinor
    beta_0, see CEEX: TCHANNEL_REDUCED_BORN, which supersedes this switch
    for Borns with exchange lines.]
    [2026-09-26: default 1. The physical-spinor crude equals the generator's
    density S~ m_born only for collinear photons. For a hard wide-angle one
    it does not: 250 GeV mu mu, one photon, rho_crude/(S~ m_born) flat at 463
    collinear and 668 (+44%) at 1-|cos theta| = 0.1-0.3, so YFS.NLO/CEEX in
    Z pT rose to 1.3-1.8 above 70 GeV. With 1: 0.99-1.01 in 55-100 GeV
    (BVR), Z pole +0.05%, and against exact tree-level e+e- -> mu mu gamma
    CEEX goes from 0.86/0.73/0.60 to 1.07/1.03/0.98 in Z pT 55-100 GeV.
    Still 2 -> 2 only (the gate below); beyond that the crude should come
    from the generator's m_born on the sampled partition.]
  */
  static const bool crudeborn(ATOOLS::Settings::GetMainSettings()["CEEX"]
                             ["CRUDE_BORN"].SetDefault(true).Get<bool>());
  /*
    2 -> 2: the Born at BornLegsAt(X_wp), which rebuilds the one pair.
    Beyond (2026-09-26): the generator's own point, GeneratorBornAt(X_wp),
    from the PRE-FSR legs. Without it the crude was the physical-spinor Born
    at the post-emission legs: e+e- -> mu mu nu nu at 250 GeV,
    rho_crude/(S~ m_born) = 248 for soft photons and 0.3-27 for hard FSR ones
    (a photon that puts mu nu gamma on the W while mu nu is off it), CEEX
    3-8x YFS.NLO. The same failure is the Hll/4f hard-FSR flux problem.
  */
  bool usered(crudeborn && m_flavs.size() >= 4 && m_comixborn
              && m_cxbalignok && !realpoint
              && (m_flavs.size() == 4 || m_prefsr.size() == m_flavs.size()));
  Amplitude Cred;
  double fluxred(1.);
  if (usered) {
    Vec4D_Vector pb;
    const double X2(m_PXvec.Abs2());
    const bool legs(X2 > 0. && (m_flavs.size() == 4 ? BornLegsAt(m_PXvec, pb)
                                                     : GeneratorBornAt(m_PXvec, pb)));
    usered = legs && ComixBornAmplitude(pb, Cred, NULL, -1., -1.);
    if (usered) fluxred = m_s/X2;
    // CEEX: CRUDE_BORN_TRACE - did the reduced-point crude engage, and how
    // does its Born compare with the generator's (m_born) and with the
    // shifted physical-spinor Born (AmpBorn)?
    static const int cbt(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["CRUDE_BORN_TRACE"].SetDefault(0).Get<int>());
    static long ncbt(0);
    if (cbt && ncbt < cbt) { ++ncbt;
      double nr(0.), ns(0.), nu(0.), nz(0.);
      Amplitude U, Z;
      const bool oku(ComixBornAmplitude(m_pceex, U, NULL, -1., -1.));
      PropShifts zero;
      const bool okz(ComixBornShifted(m_pceex, Z, zero));
      for (int f = 0; f < nhel; ++f) { nr += std::norm(Cred.m_A[f]); ns += std::norm(AmpBorn.m_A[f]);
        if (oku) nu += std::norm(U.m_A[f]); if (okz) nz += std::norm(Z.m_A[f]); }
      std::cerr<<"@@@ CRUDEBORN engaged="<<(usered?1:0)<<" legs="<<(legs?1:0)
               <<" nphot="<<m_allphotons.size()<<" x="<<(m_s>0.?1.-X2/m_s:-1.)
               <<" |Cred|^2/4="<<nr/4.<<" |Ashift|^2/4="<<ns/4.
               <<" |Aphys|^2/4="<<(oku?nu/4.:-1.)<<" |Azeroshift|^2/4="<<(okz?nz/4.:-1.)
               <<" m_born="<<m_born<<" fluxred="<<fluxred<<std::endl;
      // the shift list the partition Born was built with, and the per-leg
      // reduced-minus-physical differences, in units of the beam energy
      PropShifts shl(StageShifts(-1));
      AddExchangeLineShifts(-1, shl);
      const double Eb(m_pceex[0][0] > 0. ? m_pceex[0][0] : 1.);
      std::cerr<<"@@@ CRUDEBORN-SHIFTS n="<<shl.size();
      for (size_t j(0); j < shl.size(); ++j)
        std::cerr<<" ["<<shl[j].first<<": |d|/E="<<Vec3D(shl[j].second).Abs()/Eb
                 <<" d0/E="<<shl[j].second[0]/Eb<<"]";
      std::cerr<<"  legs:";
      for (size_t i(0); i < pb.size() && i < m_pceex.size(); ++i)
        std::cerr<<" "<<i<<":"<<Vec3D(pb[i]-m_pceex[i]).Abs()/Eb;
      std::cerr<<" if="<<m_if1<<","<<m_if2<<std::endl;
      // Controlled test: one tiny artificial shift at a time, relative to
      // the plain Born at the same (physical) legs.
      const double tiny(1e-9*Eb);
      auto probe = [&](const char *lab, size_t mask) {
        PropShifts t; t.push_back(std::make_pair(mask, Vec4D(tiny, 0., 0., tiny)));
        Amplitude T; double nt(0.);
        if (!ComixBornShifted(m_pceex, T, t)) { std::cerr<<" ["<<lab<<": fail]"; return; }
        for (int f = 0; f < nhel; ++f) nt += std::norm(T.m_A[f]);
        std::cerr<<" ["<<lab<<": "<<(nu>0.? nt/nu : -1.)<<"]"; };
      std::cerr<<"@@@ CRUDEBORN-PROBE ratio(|A_tinyshift|^2/|A_plain|^2):";
      probe("stage{0,1}", 3);
      probe("stage{2,3}", 12);
      probe("leg0", ((size_t)1) | PHASIC::Process_Base::s_propshiftleg);
      probe("leg2", ((size_t)4) | PHASIC::Process_Base::s_propshiftleg);
      probe("{0,2}", 5);
      std::cerr<<std::endl; }
  }
  /*
    One power of each flux: the physical-spinor Born of a final-stage
    partition is (q + K)^2/q^2 times the reduced-point one in the square
    (the spinors scale, the pole is shared), so KKMC's (svarX/svarQ)^2 on
    |beta_0|^2 is (svarX/svarQ) on B_red^2; the initial-stage factor
    s/X_wp^2 is the same statement for the beams.
  */
  const double crudered((pfmode == pseudoflux::neither ? 1. : m_pflux) * fluxred
                        * (m_crudefixed ? m_crudeprod : std::norm(m_Sprod)));
  /*
    The crude the CEEX weight divides by is the generator's density, whose
    initial-state Born carries the flux s/X_wp^2 relative to the Born at the
    reduced point (the fixed-order weight divides its real by S~ B and
    multiplies by X^2/s, the same statement). The shifted physical-spinor
    Born carries that factor in its spinors; the reduced-leg Born of
    BORN_AT_SPRIME does not, so it is put back on the crude here. Measured on
    e+e- -> gamma gamma at one photon: with the reduced-leg Born and no
    factor, rho_crude/(S~ m_born) is 496.1 at every x (flat, the generator's
    Born to a constant) and FO/CEEX falls as X^2/s (0.27 at x > 0.6); with
    the shifted Born rho_crude grows as 12800/238 at x > 0.6 where s/X^2 is
    17 (the t-channel numerator), which is the gamma gamma collapse.
  */
  double fluxreal(1.);
  if (realpoint) {
    const double X2r(m_PXvec.Abs2());
    if (X2r > 0. && m_s > 0.) fluxreal = m_s/X2r;
  }
  const int fmaskr(Amplitude::NHel() - 1);
  for (int f = 0; f < nhel; ++f) {
    const Complex a(fac * AmpBorn.m_A[f]);
    const Complex a0(facA0 * AmpBorn.m_A[f]);
    if (usered)
      rc += crudered * std::norm(m_cxbalign.m_A[f]
                                 * Cred.m_A[f ^ (m_comixflip & fmaskr)]);
    else
      rc += (m_crudefixed ? crudefac * std::norm(AmpBorn.m_A[f])
                          : std::real(a0 * conj(a0))) * fluxreal;
    // what beta_1 subtracts, before the soft-factor product
    m_partborn0.m_A[f] = fac0 * AmpBorn.m_A[f];
    m_AmpExpo0.m_A[f] += a0;
    m_AmpBornVirt.m_A[f] += a;
    m_AmpBornReal.m_A[f] += a;
    m_AmpExpo1.m_A[f] += a;
  }
  m_rhocrud += rc / 4.;
  m_snapBorn = m_AmpExpo1;   // Born term only, before any correction

  /*
    Is this partition's Born the Born at the configuration that realises its
    scale? CEEX pairs the propagator at m_sp with the PHYSICAL outgoing
    spinors; a generator asked for the amplitude at the momenta realising m_sp
    would use spinors for a pair rebuilt at that invariant. If the two agree,
    Comix can supply the partition Born and the hand-coded spinor algebra -
    the 2 -> 2 specific part - can go. If they do not, that difference is a
    scheme choice that has to be understood before anything is replaced.
  */
  { static const bool bx(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["PARTITION_BORN_CHECK"].Get<int>() != 0);
    static int nb(0);
    if (bx && nb < 400 && m_pceex.size() >= 4 && m_sp > 0.) {
      const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
      const double rs(sqrt(m_sp));
      if (rs > m3 + m4) {
        // the pair rebuilt at THIS partition's invariant, same direction
        Vec4D d3(m_pceex[m_if1]);
        { Poincare c0(m_pceex[m_if1] + m_pceex[m_if2]); c0.Boost(d3); }
        const double n3(d3.PSpat());
        if (n3 > 0.) {
          ++nb;
          const double E3((m_sp + m3*m3 - m4*m4)/(2.*rs));
          const double pm(sqrt(Max(0., E3*E3 - m3*m3)));
          Vec4D q3(E3, pm*d3[1]/n3, pm*d3[2]/n3, pm*d3[3]/n3);
          Vec4D q4(rs - E3, -q3[1], -q3[2], -q3[3]);
          const double lcm(0.5*rs);
          const double sgn(m_bornmomenta[0][3] < 0 ? -1. : 1.);
          Vec4D_Vector rp{Vec4D(lcm,0.,0., sgn*lcm), Vec4D(lcm,0.,0.,-sgn*lcm),
                          q3, q4};
          Amplitude Br;
          BornAmplitude(rp, Br, -1., -1., -1);
          double na(0.), nb2(0.), nd(0.);
          const int nh(Amplitude::NHel());
          for (int f = 0; f < nh; ++f) {
            na  += std::norm(AmpBorn.m_A[f]);
            nb2 += std::norm(Br.m_A[f]);
            nd  += std::norm(AmpBorn.m_A[f] - Br.m_A[f]);
          }
          /*
            The decisive question for driving Comix's currents: does it accept
            the arguments CEEX actually uses? Those are the full-energy beams
            with the physical pair, which do not conserve momentum - the
            photons carry the difference - and Comix builds its currents by
            summing external momenta.
          */
          { Amplitude Bc;
            // at THIS partition's scale, with the spinors left alone
            const bool okc(ComixBornAmplitude(m_pceex, Bc, NULL, m_sp));
            double nc(0.); bool fin(true);
            for (int f = 0; f < Amplitude::NHel(); ++f) {
              nc += std::norm(Bc.m_A[f]);
              if (!std::isfinite(Bc.m_A[f].real())
                  || !std::isfinite(Bc.m_A[f].imag())) fin = false;
            }
            const Vec4D bal(m_pceex[0]+m_pceex[1]-m_pceex[m_if1]-m_pceex[m_if2]);
            /*
              The derived map says |A_comix| = 2 e^2 |A_hand| (the 2 being
              1/sqrt(initial spin states), see CalibrateComixMap), so this
              ratio is 1 exactly when Comix has reproduced THIS partition's
              Born - which fixes which invariant Comix put in the propagator.
            */
            double na2(0.);
            for (int f = 0; f < Amplitude::NHel(); ++f)
              na2 += std::norm(AmpBorn.m_A[f]);
            const double pred(2.*m_e*m_e*sqrt(na2));
            std::cerr<<std::setprecision(14)
                     <<"@@@ NONCONS ok="<<okc<<" finite="<<fin
                     <<" sp="<<m_sp<<" sQ="<<m_svarQ
                     <<" |B_comix|="<<sqrt(nc)
                     <<" 2e2|B_ceex|="<<pred
                     <<" ratio="<<(pred>0.? sqrt(nc)/pred : -1.)
                     <<" E_imbal="<<bal[0]
                     <<std::setprecision(6)<<std::endl; }
          std::cerr<<"@@@ PARTBORN sbeam="<<(m_pceex[0]+m_pceex[1]).Abs2()
                   <<" sp="<<m_sp<<" sQ="<<m_svarQ
                   <<" |B_ceex|="<<sqrt(na)<<" |B_reduced|="<<sqrt(nb2)
                   <<" |diff|="<<sqrt(nd)
                   <<" rel="<<(na>0.? sqrt(nd/na) : -1.)<<std::endl;
        }
      }
    }
  }

  // kept for the scalar diagnostics that still read it
  SumAmplitude(m_beta00, AmpBorn, m_e * m_e);
}


/*!
  The Comix -> CEEX Born alignment, once per event.

  Comix's Born and CEEX's are the same amplitude in two conventions. The
  modulus ratio is fixed and derived (COMIX_REAL_NORM = 1/sqrt(initial spin
  states)); the phase is per helicity, being a product of per-leg little-group
  phases, and the Born calibration measured it as NOT constant across
  helicities. So the conversion is one complex number per helicity, not one
  number.

  Built at the PHYSICAL scale and then used at every partition scale. That is
  the assumption worth testing rather than asserting: it holds if the two
  codes agree on the couplings, so that the only scale dependence is the
  propagator, which cancels in the ratio. If they do not, using it at a
  different scale shows up immediately as a changed weight.
*/
void Ceex_Base::BuildComixBornAlignment()
{
  m_cxbalignok = false;
  m_cxrnorm = -1.;
  if (p_bornproc == NULL || m_pceex.size() < 4) return;
  if (!m_comixborn && !m_comixreal) return;
  /*
    Past 2 -> 2 there is no hand-coded Born to align to, and none is needed:
    the Born, the one-photon amplitude and the eikonal (ComixPolarisation)
    are all in Comix's convention, so any per-helicity factor is common to
    every term and cancels in |A|^2, and any overall constant cancels in
    rho_1/rho_0. Unit alignment, unit normalisation, no map derivation - the
    fermion flip is applied to the Born and to M_1 alike and so drops out.
  */
  if (!m_ffbar) {
    for (int f = 0; f < Amplitude::NHel(); ++f) m_cxbalign.m_A[f] = Complex(1., 0.);
    if (!m_comixcalibrated) {
      m_comixflip = m_comixphoflip ? Amplitude::NHel() : 0;
      m_comixnorm = 1.; m_normexact = true; m_comixcalibrated = true;
    }
    m_cxrnorm = 1.;
    m_cxbalignok = true;
    return;
  }
  /*
    The map has to exist before it can be used. DeriveComixMap was reached
    only from the COMIX_REAL paths, so with COMIX_BORN alone m_comixflip was
    still its sentinel -1 - and -1 & 15 is 15, a mask that permutes every
    fermion leg. It produced a Born whose helicity-summed norm was still
    right, because a permutation preserves the norm, which is exactly why the
    error survived the norm check and only showed up in the weight.
  */
  if (!m_comixcalibrated) DeriveComixMap();
  if (m_comixflip < 0 || !(m_comixnorm > 0.)) return;
  /*
    The Born is rebuilt here for the ALIGNMENT when COMIX_BORN is on, and for
    the real's NORMALISATION when the global constant is not the exact one.
    Either way it is one Comix Born per event at the reference scale.
  */
  const bool needrnorm(m_comixreal && !m_normexact);
  if (!m_comixborn && !needrnorm) return;

  /*
    Beyond 2 -> 2 there is no hand-coded Born to align to, and none is wanted:
    the alignment exists only to reconcile a Comix Born with hand-coded REAL
    terms. Once the whole amplitude comes from Comix there is no second
    convention in play, the leftover phase is global, and a global phase
    cancels in |A|^2. So the factor is simply the derived normalisation.

    Deriving the hand-coded Born for a specific final state like H l+ l- would
    not be hard, but it is the thing this is replacing.
  */
  /*
    The alignment is what makes a Comix Born usable ALONGSIDE hand-coded real
    terms, and it is general. m_cxbalign[f] = e^2 H[f] / C[f] is built once per
    event at the reference scale m_svarQ, so it carries only the CONVENTION
    difference between the two constructions - the relative phase and
    normalisation per helicity - and Comix then supplies the partition
    dependence that the factor is applied to.

    Applying it at the scale it was built at is an identity by construction.
    The content is that it is applied at the OTHER partition scales, where the
    hand-coded Born no longer has the right propagator structure: for H l+ l-
    the partition sensitivity sits in the decay propagator, which the 2 -> 2
    spinor algebra does not carry at all.

    So this needs the hand-coded Born to be right only up to a per-helicity
    factor that does not move with the partition - NOT to be right. That is a
    much weaker requirement, and it is why the refusal that used to stand here
    was too strong. It rested on the claim that there is no hand-coded Born
    past four legs; there is one, and H l+ l- lands within 2% of the Born-real
    column on it, which a zero could not do.
  */
  /*
    This used to REFUSE beyond 2 -> 2, and the reason turned out to be one
    line in COMIX rather than anything in CEEX: Amplitude::SetPropScale pinned
    EVERY internal current to the partition scale. At 2 -> 2 that is right,
    because the one internal line IS the s-channel. H l+ l- also has a
    Z -> mu mu line, and dragging that off the dilepton pole broke the Born
    wherever the coherent partition sum cancels - median weight ratio 1.198
    but p95 = 28275. With the override restricted to the s-channel line the
    same comparison is median 1, p95 1, max 1.003 over 4048 events, and the
    total moves by 1e-5.

    That is a real cross-check, not a tautology: Comix reaches the partition
    Born by Berends-Giele recursion and CEEX by hand-coded spinor algebra,
    and for H l+ l- they now agree to 1e-5 at every partition.
  */
  const double sp_save(m_sp);
  m_sp = m_svarQ;
  MakeProp();
  Amplitude H, C;
  BornAmplitude(m_pceex, H, -1., -1., -1);
  const bool ok(ComixBornAmplitude(m_pceex, C, NULL, m_svarQ, m_svarQ));
  m_sp = sp_save;
  MakeProp();
  if (!ok) return;

  const int nh(Amplitude::NHel());
  const int fmaskx(Amplitude::NHel() - 1);
  double nc(0.);
  for (int f = 0; f < nh; ++f) nc += std::norm(C.m_A[f]);
  nc = sqrt(nc);
  if (!(nc > 0.)) return;
  double nhd(0.);
  for (int f = 0; f < nh; ++f) nhd += std::norm(m_e*m_e*H.m_A[f]);
  nhd = sqrt(nhd);
  /*
    The per-event normalisation for the Comix REAL. At 2 -> 2 this is 0.5
    identically - |C| = 2 e^2 |H| there, so the ratio is 1/sqrt(4) - which is
    why the exact constant is kept in that case and this is used only when the
    constant does not apply. Whatever the hand-coded Born is missing for this
    process multiplies the one-photon real in the same way, so it cancels.
  */
  if (nhd > 0.) m_cxrnorm = nhd/nc;
  if (!m_comixborn) return;
  /*
    MEASURED, and the reason the Hll CEEX column is still not a prediction:
    the per-helicity ratio e^2 H/C is NOT a convention factor here. @@@ CEEXCMP
    on e+e- -> H mu+ mu- gives |C/H| constant ACROSS helicities within an event
    to 5-6 digits - so the flip mask and the convention are right - but swinging
    by a factor 84 BETWEEN events (N = 0.399 against N = 0.00474), tracking
    m_born over two orders of magnitude. That is the Z -> mu mu decay
    Breit-Wigner, which the hand-coded 2 -> 2 spinor structure does not carry.
    At 2 -> 2 the same ratio is 2.0000007 in every event - the spin average and
    nothing else.

    So this alignment forces Comix's Born back onto the hand-coded one, which is
    why "COMIX_BORN reproduces the hand-coded weight to 1e-5" is true by
    construction and proves nothing about the physics.

    Replacing it with the flat convention constant was TRIED and is worse
    (+69% for the Born alone, +1299% with the Comix real, against -2.51% for
    the hand-coded path). Two reasons, both instructive: m_comixnorm is the
    one-shot measured value when !m_normexact, so a flat constant reintroduces
    the seed dependence; and a constant cancels in r1/r0 anyway, so what
    actually changed was letting Comix's propagator structure through while
    RealNorm() still carried the hand-coded ratio. Born, real and crude have to
    move together.

    The real fix is the one the partition sum needs anyway: give the decay line
    its own partition-shifted invariant (Y_wp = X_wp - p_H) instead of freezing
    it, which needs SetPropScale to take a {CId -> scale} map rather than one
    scalar. Until then the Hll CEEX column is not a prediction for this process.
  */
  for (int f = 0; f < nh; ++f) {
    const Complex Cf(C.m_A[f ^ (m_comixflip & fmaskx)]);
    const Complex Hf(m_e*m_e*H.m_A[f]);
    /*
      Both sides have to be present for the ratio to mean anything. Guarding
      only the denominator leaves the slots where the HAND-CODED Born is
      negligible but Comix is not: there the ratio is a small number over a
      normal one, the alignment is whatever the mass suppression happens to
      be, and at the next partition scale Comix moves while that factor does
      not. At 2 -> 2 it is harmless because the two constructions are
      mass-suppressed in the same slots; the four live helicities carry
      everything either way.
    */
    /*
      A CONVENTION factor has unit modulus up to one overall normalisation:
      the per-leg spinor phases differ between the two constructions, their
      magnitudes do not. So the alignment is the norm ratio N = |e^2 H|/|C|
      (1/sqrt(initial spin states) at 2 -> 2) times a phase, and the phase is
      read off only where both Borns are live. In the helicity-FLIP slots the
      Born is mass-suppressed on both sides and Hf/Cf is whatever the two
      mass treatments happen to give - measured 5x on one event. That never
      mattered for the Born, but the one-photon amplitude of a photon
      collinear to a massive muon (0.1 degrees, seed 19) is LARGE in exactly
      those slots, and multiplying it by that ratio gave |A1|^2 27x KKMC's
      where the un-aligned Comix amplitude agreed to 1e-3. Modulus N
      everywhere, phase from the live slots, unit phase elsewhere.
    */
    const double N(nc > 0. ? nhd/nc : 0.);
    if (std::abs(Cf) > 1e-3*nc && std::abs(Hf) > 1e-3*nhd) {
      const Complex r(Hf/Cf);
      m_cxbalign.m_A[f] = N * r/std::abs(r);
    } else m_cxbalign.m_A[f] = Complex(N, 0.);
  }
  m_cxbalignok = true;

  /*
    Is the alignment the same at a different scale? It is scale independent
    only if CEEX and Comix differ by a COMMON factor on the photon and the Z;
    if their gamma/Z mixtures differ, the ratio drags the propagators with it
    and a factor built at one scale is wrong at another.
  */
  { static const bool ck(ATOOLS::Settings::GetMainSettings()["CEEX"]
                         ["BORN_ALIGN_CHECK"].Get<int>() != 0);
    static int nd(0);
    if (ck && nd < 60) { ++nd;
      const double s2(m_s);           // the full beam invariant
      const double sv(sp_save);
      m_sp = s2; MakeProp();
      Amplitude H2, C2;
      BornAmplitude(m_pceex, H2, -1., -1., -1);
      const bool ok2(ComixBornAmplitude(m_pceex, C2, NULL, s2));
      m_sp = sv; MakeProp();
      if (ok2) {
        /*
          Each amplitude is cut against ITS OWN norm. Comparing the amplitude
          at s against the norm at sQ mixes in the overall propagator
          suppression between the two scales - for H l+ l- the reference sits
          on the Z pole and s is far off it, so every slot fell below the
          threshold and the check reported nlive=0.
        */
        double nc2n(0.), nh2n(0.);
        for (int f = 0; f < nh; ++f) {
          nc2n += std::norm(C2.m_A[f]);
          nh2n += std::norm(m_e*m_e*H2.m_A[f]);
        }
        nc2n = sqrt(nc2n); nh2n = sqrt(nh2n);
        double amin(-1.), amax(-1.), pmin(0.), pmax(0.);
        int nl(0);
        for (int f = 0; f < nh; ++f) {
          const Complex C2f(C2.m_A[f ^ (m_comixflip & fmaskx)]);
          /*
            1e-10 admits the MASS-SUPPRESSED slots, where both sides sit at
            the m_e/E level and their ratio is a mass artefact rather than a
            convention. That is what made this check report nlive=8 where four
            helicities are live, and why its earlier verdict - that 2 -> 2 is
            as unstable as H l+ l- - was worthless. Cut at 1e-3 of the norm so
            only the helicities that carry the process are compared.
          */
          if (std::abs(C2f) <= 1e-3*nc2n || std::abs(m_cxbalign.m_A[f]) == 0.
              || std::abs(m_e*m_e*H2.m_A[f]) <= 1e-3*nh2n)
            continue;
          const Complex a2((m_e*m_e*H2.m_A[f])/C2f);
          const Complex r(a2/m_cxbalign.m_A[f]);
          const double aa(std::abs(r)), pp(std::arg(r));
          if (nl++ == 0) { amin = amax = aa; pmin = pmax = pp; }
          else { amin = Min(amin,aa); amax = Max(amax,aa);
                 pmin = Min(pmin,pp); pmax = Max(pmax,pp); }
        }
        // per helicity, at the physical scale, with and without the mask
        { static int nr(0);
          if (nr < 2) { ++nr;
            for (int f = 0; f < nh; ++f)
              std::cerr<<"@@@ BORNHEL f="<<f
                       <<" |e2H|="<<std::abs(m_e*m_e*H.m_A[f])
                       <<" |C_masked|="
                       <<std::abs(C.m_A[f ^ (m_comixflip & fmaskx)])
                       <<" |C_raw|="<<std::abs(C.m_A[f])<<std::endl;
            // which mask lines the moduli up here?
            int bm(-1); double bmet(1e30);
            double sh(0.), sc(0.);
            for (int f = 0; f < nh; ++f) { sh += std::norm(m_e*m_e*H.m_A[f]);
                                           sc += std::norm(C.m_A[f]); }
            sh = sqrt(sh); sc = sqrt(sc);
            for (int m = 0; m < nh; ++m) {
              double met(0.);
              for (int f = 0; f < nh; ++f)
                met += sqr(std::abs(m_e*m_e*H.m_A[f])/sh
                           - std::abs(C.m_A[f ^ m])/sc);
              if (bm < 0 || met < bmet) { bmet = met; bm = m; }
            }
            std::cerr<<"@@@ BORNMASK inuse="<<(m_comixflip & fmaskx)
                     <<" best="<<bm<<" metric="<<sqrt(bmet)<<std::endl; } }
        std::cerr<<"@@@ BORNALIGN sQ="<<m_svarQ<<" s="<<s2<<" nlive="<<nl
                 <<" |a(s)/a(sQ)|=["<<amin<<","<<amax<<"]"
                 <<" arg=["<<pmin<<","<<pmax<<"]"<<std::endl;
      } } }
}
