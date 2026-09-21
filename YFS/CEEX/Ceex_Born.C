/*!
  \file Ceex_Born.C

  The Born spin amplitude - the spinor structures contracted with the
  coupling x propagator factors - and beta_0^0, which is that Born dressed
  with the event's S-factor product.

*/

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
            m_bornAmp.m_A[h0][h1][h2][h3] = (CouplingZ(hel1, 1) + CouplingG()) * m_T + (CouplingZ(hel1, 1) + CouplingG()) * m_U;
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
  const double m3(Mf3 >= 0. ? Mf3 : m_flavs[2].Mass());
  const double m4(Mf4 >= 0. ? Mf4 : m_flavs[3].Mass());
  static const bool exactisrmass(
      Settings::GetMainSettings()["CEEX"]["EXACT_ISR_SPINOR_MASS"]
      .SetDefault(0).Get<int>() != 0);
  // KKMC neglects the beam mass in the spinors (its Fleps convention). For
  // Bhabha that is not available: the beam and the final-state fermion are
  // the same particle, and the t-channel ties them together on one line, so
  // a massless beam next to an exactly massive final state is inconsistent.
  // The leading helicity configurations are unchanged either way - the
  // difference is O(m_e^2/s) - but the mass-suppressed ones are not.
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
         
            m_T  = T(k[2], k[0], hel3, hel1, +1, +1, m3, mi)
                 * Tp(k[1], k[3], hel2, hel4, -1, -1, mi, m4);
            m_U  = Up(k[2], k[1], hel3, hel2, +1, +1, m3, mi)
                 * U(k[0], k[3], hel1, hel4, -1, -1, mi, m4);
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
  // KKMC adds its t-channel W (GPS_BornWPlus) as a prefactor on the s-channel
  // structures, but that does NOT carry over here. The W is purely left-handed
  // and so forces hel1 = -hel2, which is the s-channel gate; a vector gamma/Z
  // does not, and the two configurations it adds outside that gate are the pure
  // t-channel ones carrying s^2/t^2 - 86% of |M|^2 at a wide angle. They need
  // their own spinor structures, which is why m_Tampt/m_Uampt exist.
  for (int j = 0; j <= 1; j++) m_TCt[j] = m_UCt[j] = Complex(0., 0.);
  if (m_bhabha) {
    for (int j = 0; j <= 1; j++) {
      const double h(1. - 2.*j);
      // Fermi statistics puts the exchange diagram in with the opposite sign.
      // The coupling assignment mirrors the s-channel (Tt with mode 1, Ut with
      // mode 0); determined against Comix, which pins it to 0.06% across
      // angle, where a swapped assignment drifts 1.4%.
      m_TCt[j] = -(CouplingZ(h, 1) * m_propZt + CouplingG() * m_propGt);
      m_UCt[j] = -(CouplingZ(h, 0) * m_propZt + CouplingG() * m_propGt);
    }
  }
  for (int h0 = 0; h0 <= 1; h0++) {
    for (int h1 = 0; h1 <= 1; h1++) {
      for (int h2 = 0; h2 <= 1; h2++) {
        for (int h3 = 0; h3 <= 1; h3++) {
          M.m_A[h0][h1][h2][h3] = m_TC[h0] * m_Tamp[h0][h1][h2][h3]
                                + m_UC[h0] * m_Uamp[h0][h1][h2][h3];
          if (m_bhabha)
            M.m_A[h0][h1][h2][h3] += m_TCt[h0] * m_Tampt[h0][h1][h2][h3]
                                   + m_UCt[h0] * m_Uampt[h0][h1][h2][h3];
          m_bornAmp.m_A[h0][h1][h2][h3] = M.m_A[h0][h1][h2][h3];
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



bool Ceex_Base::ComixBornAmplitude(const Vec4D_Vector &p, Amplitude &A)
{
  if (p_bornproc == nullptr) return false;
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  const size_t nin(p_bornproc->NIn());
  ampl->SetNIn(nin);
  ampl->SetMS(p_bornproc->Generator());
  ampl->SetMuF2(sqr(rpa->gen.Ecms()));
  ampl->SetMuR2(sqr(rpa->gen.Ecms()));
  ampl->SetMuQ2(sqr(rpa->gen.Ecms()));
  ampl->SetMu2(sqr(rpa->gen.Ecms()));
  for (size_t i(0); i < p.size() && i < p_bornproc->Flavours().size(); ++i)
    ampl->CreateLeg(i < nin ? -p[i] : p[i], p_bornproc->Flavours()[i]);
  ampl->SetProc(p_bornproc);

  ampl->Delete();
  const double diff(0.);
  std::vector<METOOLS::Spin_Amplitudes> amps;
  if (amps.empty()) return false;

  const METOOLS::Spin_Amplitudes &sa(amps[0]);
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          const size_t idx(h0 + 2*h1 + 4*h2 + 8*h3);
          A.m_A[h0][h1][h2][h3] = (idx < sa.size()) ? sa[idx] : Complex(0.,0.);
        }

  if (m_checkxs) {
    double ss(0.);
    for (size_t i(0); i < sa.size(); ++i) ss += std::norm(sa[i]);
    std::cerr<<"@@@ CEEXCOMIX sumsq="<<ss<<" diff="<<diff
             <<" ratio="<<(ss!=0.? diff/ss : 0.)
             <<" nhel="<<sa.size()<<std::endl;
  }
  return true;
}


void Ceex_Base::InfraredSubtractedME_0_0() {
  // This partition's Born, squared on its own and added to the INCOHERENT sum
  // before it goes into the coherent m_AmpExpo0. KKMC's DistCru/CrudSum.
  double rc(0.);
  Amplitude AmpBorn;
  BornAmplitude(m_pceex, AmpBorn, -1., -1., 0);
  const Complex fac(m_e * m_e * m_cfac);
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          const Complex a(fac * AmpBorn.m_A[j1][j2][j3][j4]);
          rc += std::real(a * conj(a));
          m_AmpExpo0.m_A[j1][j2][j3][j4] += a;
          m_AmpBornVirt.m_A[j1][j2][j3][j4] += a;
          m_AmpBornReal.m_A[j1][j2][j3][j4] += a;
          m_AmpExpo1.m_A[j1][j2][j3][j4] += a;
        }
  m_rhocrud += rc / 4.;
  m_snapBorn = m_AmpExpo1;   // Born term only, before any correction

  // kept for the scalar diagnostics that still read it
  SumAmplitude(m_beta00, AmpBorn, m_e * m_e);
}
