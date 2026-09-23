/*!
  \file Ceex_Real.C

  beta_1^0: one real photon, on the initial state (KKMC HiniPlus) or the
  final state (HfinPlus), and the contractions that fold the emission
  matrices into the Born.

  Split out of Ceex_Base.C, which had reached 1867 lines. The pieces are
  grouped by what they compute, not by call order.
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


void Ceex_Base::InfraredSubtractedME_1_0(const Vec4D &k, int helk,
                                        const Complex &sProd,
                                        const Complex &Sactu, int slot, int iphot) {
  Vec4D_Vector arm1, arm2, p;
  Vec4D p1, p2, p3, p4;
  p = m_pceex;
  // helk = -1;
  // m_momenta[0],[1] are the boosted beams; m_bornmomenta is the lab copy.
  // B_{[b 1' c d]}(X): p_a replaced by the photon spinor index, and
  // B_{[1' a c d]}(X): p_b replaced. p_c,p_d are the PHYSICAL outgoing pair.
  /*
    p_c, p_d are the outgoing FERMION pair, which is legs 2 and 3 only when
    every final leg is a fermion. For e+e- -> H mu+ mu- leg 2 is the HIGGS, so
    this built the reduced Born on (gamma-index, e+, H, mu-) with mu+ absent
    entirely - and, through BornAmplitude's reduced-set fallback, with the
    HIGGS MASS in the m3 slot. The masses have to be passed explicitly for the
    same reason: a 4-vector argument makes BornAmplitude index m_flavs[2..3].
  */
  arm1 = {k, m_pceex[1], p[m_if1], p[m_if2]};
  arm2 = {m_pceex[0], k, p[m_if1], p[m_if2]};
  const double mF3(m_flavs[m_if1].Mass()), mF4(m_flavs[m_if2].Mass());
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double p1k = k * m_pceex[0];
  double p2k = k * m_pceex[1];
  Amplitude AmpBornU, AmpBornV, AmpV, AmpU;
  BornAmplitude(arm1, AmpBornU, mF3, mF4, slot >= 0 ? slot     : -1);
  BornAmplitude(arm2, AmpBornV, mF3, mF4, slot >= 0 ? slot + 1 : -1);
  const double mI(m_flavs[0].Mass());
  UGamma(k, m_pceex[0], k, helk, AmpU, 0., mI);
  VGamma(m_pceex[1], k, k, helk, AmpV, mI, 0.);

  Amplitude B10;   // per-helicity, as KKMC accumulates into m_AmpExpo1
  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) B10.m_A[Idx(a,b,c,d)] = Complex(0., 0.);

  const double gI(m_qe * m_e * m_e * m_e);
  AddU(B10, AmpBornU, AmpU,  gI / p1k / 2.);
  AddV(B10, AmpBornV, AmpV, -gI / p2k / 2.);

  m_Amp1U = AmpU;
  m_Amp1V = AmpV;

  /*
    One-shot: does the standalone reconstruction used by the per-photon Comix
    path reproduce what this routine actually computes? Compared here, in
    situ, at the same momenta, so nothing about mapping or scale can differ.
  */
  { static int nchk(0);
    if (nchk < 4 && m_perphoton) { ++nchk;
      Amplitude Hs;
      Vec4D_Vector bpc(m_pceex.begin(), m_pceex.begin()+4);
      if (HandOnePhotonAmplitude(bpc, k, helk, Hs)) {
        double nb(0.), nhs(0.), nd(0.);
        const int nhh(Amplitude::NHel());
        for (int f = 0; f < nhh; ++f) {
          nb  += std::norm(B10.m_A[f]);
          nhs += std::norm(Hs.m_A[f]);
          nd  += std::norm(Hs.m_A[f] - B10.m_A[f]);
        }
        std::cerr<<"@@@ HANDCHK |B10_isr|="<<sqrt(nb)
                 <<" |Hrecon|="<<sqrt(nhs)
                 <<" |diff|="<<sqrt(nd)
                 <<" Ek="<<k[0]<<std::endl;
      } } }

  Complex nrm(1., 0.);
  if (Sactu != Complex(0., 0.)) nrm *= sProd / Sactu;
  // Comix's amplitude in place of the hand-coded one, as a per-helicity
  // factor; see BuildComixPhotonRatios for why it is a ratio.
  const bool usecx(iphot >= 0 && iphot < (int)m_cxratiook.size()
                   && m_cxratiook[iphot]);
  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) {
          Complex v(nrm * B10.m_A[Idx(a,b,c,d)]);
          if (usecx) v *= m_cxratio[iphot].m_A[Idx(a,b,c,d)];
          m_AmpExpo1.m_A[Idx(a,b,c,d)] += v;
          m_AmpBornReal.m_A[Idx(a,b,c,d)] += v;
          if (iphot >= 0 && iphot < (int)m_realphot.size())
            m_realphot[iphot].m_A[Idx(a,b,c,d)] += v;   // <-- the result MakeRho squares
          m_snapReal.m_A[Idx(a,b,c,d)] += v;   // the beta_1^0 increment
          m_beta10 += v;                     // scalar, diagnostics only
        }
}




void Ceex_Base::InfraredSubtractedME_1_0_FSR(const Vec4D &k, int helk,
                                             const Complex &sProd,
                                             const Complex &Sactu,
                                             double CKine, int slot, int iphot) {
  const Vec4D_Vector &p(m_pceex);
  
  // The RADIATING pair, by role not by position - see the ISR routine above.
  const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
  const Vec4D_Vector arm1{p[0], p[1], k,        p[m_if2]};
  const Vec4D_Vector arm2{p[0], p[1], p[m_if1], k       };
  Amplitude AmpBornU, AmpBornV, AmpU, AmpV;
  BornAmplitude(arm1, AmpBornU, 0., m4, slot >= 0 ? slot     : -1);
  BornAmplitude(arm2, AmpBornV, m3, 0., slot >= 0 ? slot + 1 : -1);
  
  UGamma(p[m_if1], k, k, helk, AmpU, m3, 0.);
  VGamma(k, p[m_if2], k, helk, AmpV, 0., m4);

  const double p3k(k * p[m_if1]), p4k(k * p[m_if2]);
  if (IsZero(p3k) || IsZero(p4k)) return;

  Amplitude B10;
  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) B10.m_A[Idx(a,b,c,d)] = Complex(0., 0.);

  const double gF(m_qf * m_e * m_e * m_e);
  AddUF(B10, AmpU, AmpBornU,  gF / p3k / 2.);
  AddVF(B10, AmpBornV, AmpV, -gF / p4k / 2.);

  Complex nrm(1., 0.);
  if (Sactu != Complex(0., 0.)) nrm *= sProd / Sactu;

  Amplitude AmpBorn;
  BornAmplitude(m_pceex, AmpBorn, -1., -1., 0);
  const Complex kinfac(m_e * m_e * sProd * (1. - CKine));
  const bool usecxf(iphot >= 0 && iphot < (int)m_cxratiook.size()
                    && m_cxratiook[iphot]);

  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) {
          // The kinematic term is a Born piece, not part of the photon's
          // emission amplitude, so the Comix factor does not act on it.
          Complex v(nrm * B10.m_A[Idx(a,b,c,d)]);
          if (usecxf) v *= m_cxratio[iphot].m_A[Idx(a,b,c,d)];
          v += kinfac * AmpBorn.m_A[Idx(a,b,c,d)];
          m_AmpExpo1.m_A[Idx(a,b,c,d)] += v;
          m_AmpBornReal.m_A[Idx(a,b,c,d)] += v;
          if (iphot >= 0 && iphot < (int)m_realphot.size())
            m_realphot[iphot].m_A[Idx(a,b,c,d)] += v;
          m_snapReal.m_A[Idx(a,b,c,d)] += v;
          m_beta10 += v;
        }
}



void Ceex_Base::AddU(Complex &sum, const Amplitude &Born, const Amplitude &U, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * Born.m_A[Idx(j,h1,h2,h3)] * U.m_U[j][h1];
          }
          sum += CSum;
        }
      }
    }
  }
}


void Ceex_Base::AddU(Amplitude &out, const Amplitude &Born,
                     const Amplitude &U, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * Born.m_A[Idx(j,h1,h2,h3)] * U.m_U[j][h0];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}



void Ceex_Base::AddV(Complex &sum, const Amplitude &Born, const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex CSum(0, 0);
          for (int  j = 0; j <= 1; j++) {
            CSum += fac * (V.m_V[h2][j]) * Born.m_A[Idx(h0,j,h2,h3)];
          }
          sum += CSum;
        }
      }
    }
  }
}


void Ceex_Base::AddV(Amplitude &out, const Amplitude &Born,
                     const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * V.m_V[h1][j] * Born.m_A[Idx(h0,j,h2,h3)];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}

void Ceex_Base::AddUF(Amplitude &out, const Amplitude &U,
                      const Amplitude &Born, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * U.m_U[h2][j] * Born.m_A[Idx(h0,h1,j,h3)];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}


void Ceex_Base::AddVF(Amplitude &out, const Amplitude &Born,
                      const Amplitude &V, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0)
    for (int h1 = 0; h1 <= 1; ++h1)
      for (int h2 = 0; h2 <= 1; ++h2)
        for (int h3 = 0; h3 <= 1; ++h3) {
          Complex c(0., 0.);
          for (int j = 0; j <= 1; ++j)
            c += fac * Born.m_A[Idx(h0,h1,h2,j)] * V.m_V[j][h3];
          out.m_A[Idx(h0,h1,h2,h3)] += c;
        }
}




void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          sum += fac * Amp.m_A[Idx(h0,h1,h2,h3)];
        }
      }
    }
  }
}



void Ceex_Base::SumAmplitude(Complex &sum, const Amplitude &Amp1, const Amplitude &Amp2, const Complex fac) {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          for (int  j = 0; j <= 1; j++) {
            sum += fac * Amp2.m_A[Idx(j,h1,h2,h3)] * Amp1.m_U[j][h1];
          }
        }
      }
    }
  }
}
