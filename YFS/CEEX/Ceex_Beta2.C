/*!
  \file Ceex_Beta2.C

  beta_2^0 - two real photons. NOT WIRED IN: InfraredSubtractedME_2_0 is
  commented out of Calculate(), so nothing here is reachable. Kept because
  it is unfinished work rather than leftovers; see NOTES-ceex-kkmc.md.

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


Complex Ceex_Base::bsigma(Vec4D p1, Vec4D p2, int hel) {
  // eq 230 hep-ph/0006359
  Complex coeff1 = sqrt(2) * Xi(p2, p1); // sqrt(m_zeta * p2 / (m_zeta * p1));
  Vec4D phat1 = p1 - m_zeta * sqr(p1.Mass()) / (2 * m_zeta * p1);
  Vec4D phat2 = p2 - m_zeta * sqr(p2.Mass()) / (2 * m_zeta * p2);
  if (hel < 0) return coeff1 * Sminus(p1, p2);
  else return coeff1 * Splus(p1, p2);
}


void Ceex_Base::InfraredSubtractedME_2_0()
{
  int hel1, hel2, hel3, hel4;
  Complex b2double(0.0);
  Complex b2single(0.0);
  Complex b2rest(0.0);
  Complex b20(0.0);
  for (int i = 0; i < m_isrphotons.size(); ++i) {
    for (int j = 0; j < i; ++j) {
      Vec4D k1 = m_isrphotons[i];
      Vec4D k2 = m_isrphotons[j];
      b2double = BetaDouble_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2rest   = BetaRest_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2single = BetaSingle_2_0(k1, k2, m_PhoHel[i], m_PhoHel[j]);
      b2single += BetaSingle_2_0(k2, k1, m_PhoHel[j], m_PhoHel[i]);
      // b2rest = BetaRest_2_0(k1,k2,m_PhoHel[i], m_PhoHel[j]);
      b20 += b2double + b2single + b2rest;
      // b20 /= pow(2*M_PI,3)*sqrt(2);
      b20 /= pow(2 * M_PI, 6);
      b20 /= Sfactor(m_bornmomenta[0], m_bornmomenta[1], k1, m_PhoHel[i]);
      b20 /= Sfactor(m_bornmomenta[0], m_bornmomenta[1], k2, m_PhoHel[j]);

    }
  }
  m_beta20 = (b20);
}



Complex Ceex_Base::BetaDouble_2_0(Vec4D &k1, Vec4D &k2, int h1, int h2)
{
  Vec4D_Vector p;
  Vec4D p1, p2, p3, p4;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double m1 = m_momenta[0].Mass();
  double m2 = m_momenta[1].Mass();
  double m3 = m_momenta[2].Mass();
  double m4 = m_momenta[3].Mass();

  double denA = 2.*p[0] * (k1 + k2) - 2 * k1 * k2;
  double denB = 2.*p[1] * (k1 + k2) - 2 * k1 * k2;

  double deltA = 2.*k1 * k2 / denA;
  double deltB = 2.*k1 * k2 / denB;
  k1 = m_isrphotons[0];
  k2 = m_isrphotons[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], h1) / (2 * k1 * p[0]);
  Complex s2a = -m_qe * bsigma(k2, p[0], h2) / (2 * k2 * p[0]);
  Complex s1b = m_qe * bsigma(k1, p[1], h1) / (2 * k1 * p[1]);
  Complex s2b = m_qe * bsigma(k2, p[1], h2) / (2 * k2 * p[1]);
  return 4 * M_PI * m_alpha * (s1a * s2a * deltA + s1b * s2b * deltB) * BornAmplitude(p);
}



Complex Ceex_Base::BetaSingle_2_0(Vec4D &k1, Vec4D &k2, int hk1, int hk2) {
  int hel1, hel2, hel3, hel4;
  // r_if = 2k_i*p_f
  // r_ij = 2k_i*k_j
  Complex Single;
  Vec4D p1, p2, p3, p4;
  Vec4D_Vector p;
  Amplitude AmpU1, AmpU2, AmpV1, AmpV2;
  Amplitude AmpBornU1, AmpBornU2, AmpBornV1, AmpBornV2;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);
  double m1 = m_momenta[0].Mass();
  double m2 = m_momenta[1].Mass();
  double m3 = m_momenta[2].Mass();
  double m4 = m_momenta[3].Mass();
  Vec4D Q;
  Q = m_bornmomenta[0] + m_bornmomenta[1] - k1 - k2;
  double r11 = 2 * k1 * k1;
  double r12 = 2 * k1 * k2;

  double r1a = 2 * k1 * p[0];
  double r1b = 2 * k1 * p[1];

  double r2a = 2 * k2 * p[0];
  double r2b = 2 * k2 * p[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], hk1) / (2 * k1 * p[0]);
  Complex s1b = m_qe * bsigma(k1, p[1], hk1) / (2 * k1 * p[1]);

  Complex s2a = -m_qe * bsigma(k2, p[0], hk2) / (2 * k2 * p[0]);
  Complex s2b = m_qe * bsigma(k2, p[1], hk2) / (2 * k2 * p[1]);
  // BornAmplitude(arm1,AmpBornU);
  // BornAmplitude(arm2, AmpBornV);
  // UGamma(m_bornmomenta[0],k, k, helk,AmpU);
  // VGamma(k, m_bornmomenta[1], k, helk,AmpV);

  // AddU(m_beta10, AmpBornU, AmpU, m_e/p1k/2.);
  // AddV(m_beta10, AmpBornV, AmpV, -m_e/p2k/2.);
  Vec4D_Vector armU1 = {k1, m_bornmomenta[1], p[2], p[3]};
  Vec4D_Vector armU2 = {k2, m_bornmomenta[1], p[2], p[3]};
  Vec4D_Vector armV1 = {m_bornmomenta[0], k1, p[2], p[3]};
  Vec4D_Vector armV2 = {m_bornmomenta[0], k2, p[2], p[3]};

  BornAmplitude(armU2, AmpBornU2);
  UGamma(m_bornmomenta[0], k2, k1, hk2, AmpU2);
  Complex f1 = -s2a / (r12 - r1a - r2a);
  AddU(Single, AmpBornU2, AmpU2, f1);

  BornAmplitude(armV1, AmpBornV1);
  VGamma(k1, k2, m_bornmomenta[1], hk1, AmpV1);
  f1 = s1b / (r12 - r1a - r2a);
  AddV(Single, AmpBornV1, AmpV1, f1);

  BornAmplitude(armU1, AmpBornU1);
  UGamma(m_bornmomenta[0], k1, k1, hk1, AmpU1);
  f1 = -s2a * (1 / (r12 - r1a - r2a) - 1. / (-r2b));
  AddU(Single, AmpBornU1, AmpU1);



  BornAmplitude(armV2, AmpBornV2);
  VGamma(k2, k2, m_bornmomenta[1], hk1, AmpU1);
  f1 = s1b * (1 / (r12 - r1a - r2a) - 1. / (-r2b));
  AddU(Single, AmpBornU1, AmpU1);



  return Single;
}


Complex Ceex_Base::BetaRest_2_0(Vec4D &k1, Vec4D &k2, int hk1, int hk2) {
  int hel1, hel2, hel3, hel4;
  Vec4D p1, p2, p3, p4;
  Vec4D_Vector p;
  p = m_bornmomenta;
  Complex amp1ISR(0, 0), sum1(0, 0), sum2(0, 0);

  hk1 = 0.5 * (3 - hk1);
  hk2 = 0.5 * (3 - hk2);
  double r11 = 2 * k1 * k1;
  double r12 = 2 * k1 * k2;

  double r1a = 2 * k1 * p[0];
  double r1b = 2 * k1 * p[1];

  double r2a = 2 * k2 * p[0];
  double r2b = 2 * k2 * p[1];
  Complex s1a = -m_qe * bsigma(k1, p[0], hk1) / (2 * k1 * p[0]);
  Complex s2a = -m_qe * bsigma(k2, p[0], hk2) / (2 * k2 * p[0]);

  Complex s1b = m_qe * bsigma(k1, p[1], hk1) / (2 * k1 * p[1]);
  Complex s2b = m_qe * bsigma(k2, p[1], hk2) / (2 * k2 * p[1]);

  Complex betarest(0.0);
  Complex  t1, t2, t3, t4;

  Complex Amp, Single;
  Amplitude AmpU2, AmpBornU2;
  Amp = BornAmplitude(m_bornmomenta);
  UGamma(m_bornmomenta[0], k2, k1, hk2, AmpU2);
  Complex f1 = -s2a / (r12 - r1a - r2a);
  AddU(Single, AmpBornU2, AmpU2, f1);


  return betarest;///(Sfactor(m_bornmomenta[0], m_bornmomenta[1], k1, hk1)*Sfactor(m_bornmomenta[0], m_bornmomenta[1], k2, hk2));
}
