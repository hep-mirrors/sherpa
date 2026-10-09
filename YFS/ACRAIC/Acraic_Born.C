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

using namespace YFS;


Complex Acraic_Base::BornAmplitude(const Vec4D_Vector &k) {
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


// Spinor structures m_Tamp/m_Uamp are cached per call site for one event;
// only m_TC/m_UC depend on the partition. Slot 0 is the Born, photon j owns
// 1+4j..4+4j (its two ISR and two FSR arms).
void Acraic_Base::BornAmplitude(const Vec4D_Vector &k, Amplitude &M,
                             double Mf3, double Mf4, int slot) {
  Complex amp;
  int hel1, hel2, hel3, hel4;
  const double m1(m_flavs[0].Mass());
  // The outgoing fermion pair: m_if1/m_if2 in the full ACRAIC set (not 2,3:
  // for H l+ l- leg 2 is the Higgs), but always 2,3 in the reduced 2 -> 2 set
  // {P1, P2, q3, q4} the real paths pass.
  const size_t i1(k.size() > 4 ? m_if1 : 2), i2(k.size() > 4 ? m_if2 : 3);
  const double m3(Mf3 >= 0. ? Mf3 : m_flavs[i1].Mass());
  const double m4(Mf4 >= 0. ? Mf4 : m_flavs[i2].Mass());
  static const bool exactisrmass(
      Settings::GetMainSettings()["ACRAIC"]["EXACT_ISR_SPINOR_MASS"]
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
  // the only partition-dependent part: two coupling x propagator factors
  for (int j = 0; j <= 1; j++) {
    double h = 1. - 2.*j;
    m_UC[j] = CouplingZ(h, 0) * m_propZ + CouplingG() * m_propG;
    m_TC[j] = CouplingZ(h, 1) * m_propZ + CouplingG() * m_propG;
  }
  // Bhabha t-channel couplings, carried by the crossed structures built above
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







Complex Acraic_Base::BornAmplitude(const Vec4D_Vector &k, int h0, int h1, int h2, int h3) {
  Complex amp;
  if (h0 == -h1 ) {
    m_T = T(k[2], k[0], h2, h0) * Tp(k[1], k[3], h1, h3);
    m_U = Up(k[2], k[1], h2, h1) * U(k[0], k[3], h0, h3);
    amp = (CouplingZ(h0, 0) * m_propZ + CouplingG() * m_propG) * m_U + (CouplingZ(h0, 1) * m_propZ + CouplingG() * m_propG) * m_T;
  }
  return amp;
}



Complex Acraic_Base::BornAmplitude(Vec4D p1, Vec4D p2, Vec4D p3, Vec4D p4, int h0, int h1, int h2, int h3)
{
  Vec4D_Vector tmp;
  tmp.push_back(p1);
  tmp.push_back(p2);
  tmp.push_back(p3);
  tmp.push_back(p4);
  return BornAmplitude(tmp, h0, h1, h2, h3);
}



//! Born spin amplitudes from Comix; false if unavailable.
bool Acraic_Base::ComixBornAmplitude(const Vec4D_Vector &p, Amplitude &A,
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


bool Acraic_Base::ComixBornShifted(const Vec4D_Vector &p, Amplitude &A,
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

// Adds this partition's Born, e^2 (flux) beta_0(X_wp) prod_i s_i^{wp_i}, and
// stores it without the soft-factor product (m_partborn0) so beta_1 subtracts
// exactly what was added.
void Acraic_Base::InfraredSubtractedME_0_0() {
  // this partition's Born, squared on its own and added to the incoherent sum
  double rc(0.);
  Amplitude AmpBorn;
  BornAmplitude(m_pacraic, AmpBorn, -1., -1., 0);
  // Comix's Born at this partition's scale, in ACRAIC's convention: spinors
  // at the physical momenta, only the propagator pole moves.
  bool realpoint(false);
  if (m_comixborn && m_cxbalignok) {
    Amplitude C;
    /*
      beta_0(X_wp), ACRAIC: BORN_AT_SPRIME:
        0  KKMC's: physical spinors, pole pinned to m_sp (decay line to
           m_svarY), times the pseudo-flux svarY/svarQ;
        1  the Born at the real point BornLegsAt(X_wp), no flux: what Comix's
           M_1 reduces to as its photon goes soft.
    */
    static const bool sprime(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                             ["BORN_AT_SPRIME"].Get<bool>());
    bool ok(false);
    // BornLegsAt is a 2 -> 2 reduction (one pair at Y): not with W stages
    if ((sprime || m_redborn) && !WStagesActive()) {   // m_redborn: ACRAIC: TCHANNEL_REDUCED_BORN
      Vec4D_Vector pb;
      ok = BornLegsAt(m_PXvec, pb)
           && ComixBornAmplitude(pb, C, NULL, -1., -1.);
      realpoint = ok;
    }
    // The pinned Born by momentum shift: initial-state photons added to the
    // initial-side lines, final-state ones to the final side, so both give
    // X_wp^2. Same as the scalar pin (checked below) and what M_1 reduces to.
    if (!ok) {
      PropShifts sh(StageShifts(-1));
      AddExchangeLineShifts(-1, sh);     // space-like exchange lines
      ok = ComixBornShifted(m_pacraic, C, sh);
      if (ok && m_checkxs && m_ffbar && m_nstages == 2) {
        static int nchk(0);
        if (nchk < 12) { ++nchk;
          Amplitude Cp;
          if (ComixBornAmplitude(m_pacraic, Cp, NULL, m_sp, m_svarY)) {
            double nd(0.), nn(0.);
            for (int f = 0; f < Amplitude::NHel(); ++f) {
              nd += std::norm(C.m_A[f] - Cp.m_A[f]); nn += std::norm(Cp.m_A[f]); }
            std::cerr<<"@@@ SHIFTCHK nphot="<<m_allphotons.size()
                     <<" sp="<<m_sp<<" sY="<<m_svarY
                     <<" |shift-pin|/|pin|="<<(nn>0.? sqrt(nd/nn) : -1.)
                     <<std::endl; } }
      }
    }
    if (!ok) ok = ComixBornAmplitude(m_pacraic, C, NULL, m_sp, m_svarY);
    if (ok) {
      const int nh(Amplitude::NHel());
      for (int f = 0; f < nh; ++f)
        AmpBorn.m_A[f] = m_cxbalign.m_A[f]*C.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]
                         /(m_e*m_e);
    }
  }
  // The real-point Born carries its own scale, so no flux. fac0 feeds rho_1
  // and beta_1's subtraction, fac00 rho_0 (ACRAIC: NO_PSEUDOFLUX).
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  const Complex fac0(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::rho0_and_rho1 ? m_pflux : 1.)));
  const Complex fac00(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::neither ? 1. : m_pflux)));
  const Complex fac(fac0 * m_Sprod), facA0(fac00 * m_Sprod);
  // every helicity the container holds, not only the 16 of a 2 -> 2
  const int nhel(Amplitude::NHel());
  // with a collapsed photon the crude's soft product is the incoherent
  // m_crudeprod, not |m_Sprod|^2
  const double crudefac(std::norm(fac00) * m_crudeprod);
  /*
    ACRAIC: CRUDE_BORN. The crude's Born is the generator's: the Born at the
    partition's reduced point times s/X_wp^2. The coherent sum and beta_1's
    subtraction keep the physical-spinor form above, which is what M_1
    reduces to. The two agree for an s-channel Born, but not with a
    space-like exchange line nor for hard wide-angle photons.
  */
  static const bool crudeborn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                             ["CRUDE_BORN"].SetDefault(true).Get<bool>());
  // 2 -> 2: the Born at BornLegsAt(X_wp), which rebuilds the one pair;
  // beyond: the generator's point GeneratorBornAt(X_wp) from the pre-FSR legs.
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
    // ACRAIC: CRUDE_BORN_TRACE: the reduced-point crude's Born against m_born
    // and the shifted physical-spinor Born
    static const int cbt(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["CRUDE_BORN_TRACE"].SetDefault(0).Get<int>());
    static long ncbt(0);
    if (cbt && ncbt < cbt) { ++ncbt;
      double nr(0.), ns(0.), nu(0.), nz(0.);
      Amplitude U, Z;
      const bool oku(ComixBornAmplitude(m_pacraic, U, NULL, -1., -1.));
      PropShifts zero;
      const bool okz(ComixBornShifted(m_pacraic, Z, zero));
      for (int f = 0; f < nhel; ++f) { nr += std::norm(Cred.m_A[f]); ns += std::norm(AmpBorn.m_A[f]);
        if (oku) nu += std::norm(U.m_A[f]); if (okz) nz += std::norm(Z.m_A[f]); }
      std::cerr<<"@@@ CRUDEBORN engaged="<<(usered?1:0)<<" legs="<<(legs?1:0)
               <<" nphot="<<m_allphotons.size()<<" x="<<(m_s>0.?1.-X2/m_s:-1.)
               <<" |Cred|^2/4="<<nr/4.<<" |Ashift|^2/4="<<ns/4.
               <<" |Aphys|^2/4="<<(oku?nu/4.:-1.)<<" |Azeroshift|^2/4="<<(okz?nz/4.:-1.)
               <<" m_born="<<m_born<<" fluxred="<<fluxred<<std::endl;
      // the partition Born's shifts and per-leg reduced-minus-physical, / E_beam
      PropShifts shl(StageShifts(-1));
      AddExchangeLineShifts(-1, shl);
      const double Eb(m_pacraic[0][0] > 0. ? m_pacraic[0][0] : 1.);
      std::cerr<<"@@@ CRUDEBORN-SHIFTS n="<<shl.size();
      for (size_t j(0); j < shl.size(); ++j)
        std::cerr<<" ["<<shl[j].first<<": |d|/E="<<Vec3D(shl[j].second).Abs()/Eb
                 <<" d0/E="<<shl[j].second[0]/Eb<<"]";
      std::cerr<<"  legs:";
      for (size_t i(0); i < pb.size() && i < m_pacraic.size(); ++i)
        std::cerr<<" "<<i<<":"<<Vec3D(pb[i]-m_pacraic[i]).Abs()/Eb;
      std::cerr<<" if="<<m_if1<<","<<m_if2<<std::endl;
      // one tiny artificial shift at a time, relative to the plain Born
      const double tiny(1e-9*Eb);
      auto probe = [&](const char *lab, size_t mask) {
        PropShifts t; t.push_back(std::make_pair(mask, Vec4D(tiny, 0., 0., tiny)));
        Amplitude T; double nt(0.);
        if (!ComixBornShifted(m_pacraic, T, t)) { std::cerr<<" ["<<lab<<": fail]"; return; }
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
  // One power of each flux: |B_phys|^2 = (q+K)^2/q^2 |B_red|^2 for a final
  // stage, so KKMC's (svarX/svarQ)^2 on |beta_0|^2 is svarX/svarQ on B_red^2;
  // s/X_wp^2 is the same for the beams.
  const double crudered((pfmode == pseudoflux::neither ? 1. : m_pflux) * fluxred
                        * (m_crudefixed ? m_crudeprod : std::norm(m_Sprod)));
  // The generator's density carries s/X_wp^2 relative to the reduced-point
  // Born. The shifted physical-spinor Born has it in its spinors; the
  // BORN_AT_SPRIME real-point Born does not, so it is put on the crude here.
  double fluxreal(1.);
  if (realpoint) {
    const double X2r(m_PXvec.Abs2());
    if (X2r > 0. && m_s > 0.) fluxreal = m_s/X2r;
  }
  for (int f = 0; f < nhel; ++f) {
    const Complex a(fac * AmpBorn.m_A[f]);
    const Complex a0(facA0 * AmpBorn.m_A[f]);
    if (usered)
      rc += crudered * std::norm(m_cxbalign.m_A[f]
                                 * Cred.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
    else
      rc += (m_crudefixed ? crudefac * std::norm(AmpBorn.m_A[f])
                          : std::real(a0 * conj(a0))) * fluxreal;
    // what beta_1 subtracts, before the soft-factor product
    m_partborn0.m_A[f] = fac0 * AmpBorn.m_A[f];
    m_AmpExpo0.m_A[f] += a0;
    m_AmpBornVirt.m_A[f] += a;
    m_AmpBornOnly.m_A[f] += a;
    m_AmpBornReal.m_A[f] += a;
    m_AmpExpo1.m_A[f] += a;
  }
  m_rhocrud += rc / 4.;
  m_snapBorn = m_AmpExpo1;   // Born term only, before any correction

  CheckPartitionBorn(AmpBorn);

  // kept for the scalar diagnostics that still read it
  SumAmplitude(m_beta00, AmpBorn, m_e * m_e);
}


/*!
  The Comix -> ACRAIC Born alignment, once per event: one complex factor per
  helicity (per-leg little-group phases), modulus COMIX_REAL_NORM =
  1/sqrt(initial spin states) at 2 -> 2. Built at the physical scale and used
  at every partition scale. When that constant is not exact (m_normexact) the
  real uses m_cxrnorm, the same ratio rebuilt from this event's Born.
*/
void Acraic_Base::BuildComixBornAlignment()
{
  m_cxbalignok = false;
  m_cxrnorm = -1.;
  if (p_bornproc == NULL || m_pacraic.size() < 4) return;
  if (!m_comixborn && !m_comixreal) return;
  // No hand-coded Born beyond f fbar: Born, M_1 and eikonal (ComixPolarisation)
  // are all in Comix's convention, so unit alignment and normalisation.
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
  // Calibrate first: the sentinel m_comixflip = -1 would act as a mask
  // flipping every fermion leg, which keeps the norm and so passes unnoticed.
  if (!m_comixcalibrated) DeriveComixMap();
  if (m_comixflip < 0 || !(m_comixnorm > 0.)) return;
  // one Comix Born per event at the reference scale, for the alignment
  // (COMIX_BORN) and/or the real's per-event normalisation
  const bool needrnorm(m_comixreal && !m_normexact);
  if (!m_comixborn && !needrnorm) return;

  // m_cxbalign[f] = e^2 H[f]/C[f] at the reference scale m_svarQ carries only
  // the convention difference; Comix supplies the partition dependence, so H
  // need only be right up to a partition-independent factor per helicity.
  const double sp_save(m_sp);
  m_sp = m_svarQ;
  MakeProp(false);   // like the unshifted Comix Born it is compared with
  Amplitude H, C;
  BornAmplitude(m_pacraic, H, -1., -1., -1);
  const bool ok(ComixBornAmplitude(m_pacraic, C, NULL, m_svarQ, m_svarQ));
  m_sp = sp_save;
  MakeProp();
  if (!ok) return;

  const int nh(Amplitude::NHel());
  double nc(0.);
  for (int f = 0; f < nh; ++f) nc += std::norm(C.m_A[f]);
  nc = sqrt(nc);
  if (!(nc > 0.)) return;
  double nhd(0.);
  for (int f = 0; f < nh; ++f) nhd += std::norm(m_e*m_e*H.m_A[f]);
  nhd = sqrt(nhd);
  // per-event normalisation of the Comix real (0.5 identically at 2 -> 2,
  // where the exact constant is used instead)
  if (nhd > 0.) m_cxrnorm = nhd/nc;
  if (!m_comixborn) return;
  // Meaningful only while e^2 H/C is a pure convention factor; where H lacks
  // structure Comix has (e.g. a decay Breit-Wigner) this forces Comix's Born
  // back onto H.
  for (int f = 0; f < nh; ++f) {
    const Complex Cf(C.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
    const Complex Hf(m_e*m_e*H.m_A[f]);
    // Modulus N everywhere; the phase only where both Borns are live, since
    // in mass-suppressed slots Hf/Cf is a mass artefact that collinear
    // photons amplify.
    const double N(nc > 0. ? nhd/nc : 0.);
    if (std::abs(Cf) > 1e-3*nc && std::abs(Hf) > 1e-3*nhd) {
      const Complex r(Hf/Cf);
      m_cxbalign.m_A[f] = N * r/std::abs(r);
    } else m_cxbalign.m_A[f] = Complex(N, 0.);
  }
  m_cxbalignok = true;

  CheckBornAlignment(H, C, nh, sp_save);
}
