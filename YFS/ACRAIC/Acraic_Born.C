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


void Acraic_Base::BornAmplitude(const Vec4D_Vector &k, Amplitude &M,
                             double Mf3, double Mf4, int slot) {
  Complex amp;
  int hel1, hel2, hel3, hel4;
  const double m1(m_flavs[0].Mass());
  const size_t i1(k.size() > 4 ? m_if1 : 2), i2(k.size() > 4 ? m_if2 : 3);
  const double m3(Mf3 >= 0. ? Mf3 : m_flavs[i1].Mass());
  const double m4(Mf4 >= 0. ? Mf4 : m_flavs[i2].Mass());
  const double mi(m_bhabha ? m1 : 0.);
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
        }
      }
    }
  }
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
  // FF_SCALE_PROBE (debug): the shift list handed to Comix
  static const bool probe(ATOOLS::Settings::GetMainSettings()["FF_SCALE_PROBE"]
                          .SetDefault(false).Get<bool>());
  static int nprobe(0);
  if (probe && nprobe < 30 && p.size() > 3 && (p[2] + p[3]).Abs2() < 0.95*(p[0] + p[1]).Abs2()) {
    ++nprobe;
    std::cerr<<"@@@BSHIFT nphot="<<m_allphotons.size()<<" nstage="<<m_stage.size()
             <<" nshifts="<<shifts.size()<<" (pi pair)^2/s="
             <<(p.size() > 3 ? (p[2] + p[3]).Abs2()/(p[0] + p[1]).Abs2() : -1.);
    for (size_t g(0); g < shifts.size(); ++g)
      std::cerr<<"  [mask "<<shifts[g].first<<" shift "<<shifts[g].second<<"]";
    std::cerr<<std::endl;
  }
  std::vector<METOOLS::Spin_Amplitudes> amps;
  if (!p_bornproc->BornSpinAmplitudesShifts(p, amps, NULL, shifts)) return false;
  if (amps.empty()) return false;
  const METOOLS::Spin_Amplitudes &sa(amps[0]);
  const int nh(Amplitude::NHel());
  if ((int)sa.size() < nh) return false;
  for (int f = 0; f < nh; ++f) A.m_A[f] = sa[f];
  return true;
}

void Acraic_Base::InfraredSubtractedME_0_0() {
  double rc(0.);
  Amplitude AmpBorn;
  BornAmplitude(m_pacraic, AmpBorn, -1., -1., 0);
  bool realpoint(false);
  PropShifts bornsh;   // the shifts of the Comix Born, if it took them
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
      if (ok) bornsh = sh;
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
  // ACRAIC: TNP_IFI: the Born of the flipped weight, greco's a_g -> -a_g
  Amplitude AmpBornFlip(AmpBorn);
  if (m_tnpifi && !m_resexp.empty()) {
    Amplitude Cb;
    if (!bornsh.empty() && ComixBornShifted(m_pacraic, Cb, TnpIfiFlippedResonance(bornsh))) {
      for (int f = 0; f < Amplitude::NHel(); ++f)
        AmpBornFlip.m_A[f] = m_cxbalign.m_A[f]*Cb.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]
                             /(m_e*m_e);
    }
    else if (!(m_comixborn && m_cxbalignok)) ++m_tnpifievenres;
  }
  static const pseudoflux::code pfmode(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["NO_PSEUDOFLUX"].Get<pseudoflux::code>());
  const Complex fac0(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::rho0_and_rho1 ? m_pflux : 1.)));
  const Complex fac00(m_e * m_e * (realpoint ? 1. : (pfmode == pseudoflux::neither ? 1. : m_pflux)));
  const Complex fac(fac0 * m_Sprod), facA0(fac00 * m_Sprod);
  // every helicity the container holds, not only the 16 of a 2 -> 2
  const int nhel(Amplitude::NHel());
  const double crudefac(std::norm(fac00) * m_crudeprod);
  
  static const bool crudeborn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                             ["CRUDE_BORN"].SetDefault(true).Get<bool>());
  static const crudebornlegs::code cbl(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                                       ["CRUDE_BORN_LEGS"].Get<crudebornlegs::code>());
  // ACRAIC: CRUDE_BORN_LEGS stage: the partition's own channel point
  // (StageBornAt) beyond 2 -> 2; 2 -> 2 keeps BornLegsAt, which is the same point
  const bool stagelegs(cbl == crudebornlegs::stage && m_flavs.size() > 4);
  const bool genlegs(!stagelegs && (cbl == crudebornlegs::generator || m_flavs.size() > 4)
                     && m_prefsr.size() == m_flavs.size());
  bool usered(crudeborn && m_flavs.size() >= 4
              && m_cxbalignok && !realpoint
              && (m_flavs.size() == 4 || genlegs || stagelegs));
  Amplitude Cred;
  double fluxred(1.);
  if (usered) {
    Vec4D_Vector pb;
    const double X2(m_PXvec.Abs2());
    const bool legs(X2 > 0. && (stagelegs ? StageBornAt(m_PXvec, pb)
                                : genlegs ? GeneratorBornAt(m_PXvec, pb)
                                          : BornLegsAt(m_PXvec, pb)));
    usered = legs && ComixBornAmplitude(pb, Cred, NULL, -1., -1.);
    if (usered) fluxred = m_s/X2;
  }
  const double crudered((pfmode == pseudoflux::neither ? 1. : m_pflux * m_pfluxmass) * fluxred
                        * (m_crudefixed ? m_crudeprod : std::norm(m_Sprod)));
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
    TnpAdd(tnp_b0, f, a);
    // ACRAIC: TNP_IFI: the soft product with the final-stage charges flipped
    if (m_tnpifi) {
      m_partborn0flip.m_A[f] = fac0 * AmpBornFlip.m_A[f];
      TnpIfiAdd(f, m_tnppartflip*facA0*AmpBornFlip.m_A[f], m_tnppartflip*fac*AmpBornFlip.m_A[f]);
    }
  }
  m_rhocrud += rc / 4.;
  m_snapBorn = m_AmpExpo1;   // Born term only, before any correction
}


void Acraic_Base::BuildComixBornAlignment()
{
  m_cxbalignok = false;
  m_cxrnorm = -1.;
  if (p_bornproc == NULL || m_pacraic.size() < 4) return;
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
  if (!m_comixcalibrated) DeriveComixMap();
  if (m_comixflip < 0 || !(m_comixnorm > 0.)) return;
  // one Comix Born per event at the reference scale, for the alignment
  // (COMIX_BORN) and/or the real's per-event normalisation

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
  // Meaningful only while e^2 H/C is a pure convention factor; where H lacks
  // structure Comix has (e.g. a decay Breit-Wigner) this forces Comix's Born
  // back onto H.
  for (int f = 0; f < nh; ++f) {
    const Complex Cf(C.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
    const Complex Hf(m_e*m_e*H.m_A[f]);
    const double N(nc > 0. ? nhd/nc : 0.);
    if (std::abs(Cf) > 1e-3*nc && std::abs(Hf) > 1e-3*nhd) {
      const Complex r(Hf/Cf);
      m_cxbalign.m_A[f] = N * r/std::abs(r);
    } else m_cxbalign.m_A[f] = Complex(N, 0.);
  }
  m_cxbalignok = true;
}
