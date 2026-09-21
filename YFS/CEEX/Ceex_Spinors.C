/*!
  \file Ceex_Spinors.C

  The GPS spinor algebra of hep-ph/0006359: the T/U structures, the basic
  spinor products, the photon emission matrices U/V, and the soft factor.
  Nothing here knows about partitions, betas or couplings.

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



Complex Ceex_Base::T(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                     double M1, double M2) {
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = Xi(p1, p2); //sqrt(m_zeta*p1/(m_zeta*p2));
  double sq2 = Xi(p2, p1); //sqrt(m_zeta*p2/(m_zeta*p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else if (h1 == -1) s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    s = double(C1) * m1 * sq2 + double(C2) * m2 * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}



Complex Ceex_Base::Tp(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                      double M1, double M2) {
  h1 = -h1;
  h2 = -h2;
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = Xi(p1, p2); //sqrt(m_zeta*p1/(m_zeta*p2));
  double sq2 = Xi(p2, p1); //sqrt(m_zeta*p2/(m_zeta*p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    s = double(C1) * m1 * sq2 + double(C2) * m2 * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::U(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                     double M1, double M2) {
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = sqrt(m_zeta * p1 / (m_zeta * p2));
  double sq2 = sqrt(m_zeta * p2 / (m_zeta * p1));
  if (h1 == -h2) {
    if (h1 == -1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s = double(C1) * m1 * sq2 + double(C2) * m2 * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::Up(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                      double M1, double M2) {
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = sqrt(m_zeta * p1 / (m_zeta * p2));
  double sq2 = sqrt(m_zeta * p2 / (m_zeta * p1));
  if (h1 == -h2) {
    if (h1 == 1) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
    return s;
  }
  else if (h1 == h2) {
    s = double(C1) * m1 * sq2 + double(C2) * m2 * sq1;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::S(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
  Complex s(0, 0);
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s = p1.Mass() * sq1 + p2.Mass() * sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}


Complex Ceex_Base::S(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2) {
  Complex s(0, 0);
  Vec4D phat1, phat2;
  if (IsEqual(m1, 0)) h1 = -h1;
  if (IsEqual(m2, 0)) h2 = -h2;
  double sq1 = sqrt(m_zeta * p2 / (m_zeta * p1));
  double sq2 = sqrt(m_zeta * p1 / (m_zeta * p2));
  if (h1 == -h2) {
    if (h1 > 0) s = Splus(p1, p2);
    else s = Sminus(p1, p2);
  }
  else if (h1 == h2) {
    s = m1 * sq1 + m2 * sq2;
  }
  else {
    msg_Error() << METHOD << "Wrong helicities\n";
  }
  return s;
}



Complex Ceex_Base::Splus(const Vec4D &p, const Vec4D &q) {
  Complex sp = -Complex(q[2], q[3]) * sqrt((p[0] - p[1]) / (q[0] - q[1]));
  sp += Complex(p[2], p[3]) * sqrt((q[0] - q[1]) / (p[0] - p[1]));
  return sp;
}


Complex Ceex_Base::Sminus(const Vec4D &p, const Vec4D &q) {
  return -conj(Splus(p, q));
}



void Ceex_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU,
                       double M1, double M2) {

  AmpU.m_U[0][0] = UGamma(p1, p2, k, sigma,  1,  1, false, M1, M2);
  AmpU.m_U[0][1] = UGamma(p1, p2, k, sigma,  1, -1, false, M1, M2);
  AmpU.m_U[1][0] = UGamma(p1, p2, k, sigma, -1,  1, false, M1, M2);
  AmpU.m_U[1][1] = UGamma(p1, p2, k, sigma, -1, -1, false, M1, M2);

}



Complex Ceex_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j,
                          bool negMass, double M1, double M2) {
  // eq 222 hep-ph/0006359
  if (h1 == -1) return -conj(UGamma(p2, p1, k, 1, j, i, negMass, M2, M1));
  double sqr2 = sqrt(2);
  double m1 = (M1 >= 0. ? M1 : p1.Mass());
  double m2 = (M2 >= 0. ? M2 : p2.Mass());
  if (negMass) {
    m1 = -m1;
    m2 = -m2;
  }
  Vec4D p1hat = p1 - m_zeta * m1 * m1 / (2 * m_zeta * p1);
  Vec4D p2hat = p2 - m_zeta * m2 * m2 / (2 * m_zeta * p2);
  Complex coeff1 = sqrt(2) * Xi(p1, k); // sqrt(m_zeta * p1 / (m_zeta * k));
  Complex coeff2 = sqrt(2) * Xi(p2, k); // sqrt(m_zeta * p2 / (m_zeta * k));
  if ( i == 1 && j == 1) return sqr2 * Xi(p2, k) * Splus(k, p1hat);
  else if (i == 1 && j == -1) return 0;
  else if (i == -1 && j == 1) return sqr2 * (m2 * Xi(p1, p2) - m1 * Xi(p2, p1));
  else if (i == -1 && j == -1) return sqr2 * Xi(p1, k) * Splus(k, p2hat);
  else {
    msg_Out() << "h1 = " << i << "\n h2 = " << j;
    THROW(fatal_error, "Wrong helicities in CEEX");
  }
}


void Ceex_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV,
                       double M1, double M2) {

  AmpV.m_V[0][0] = VGamma(p1, p2, k, sigma,  1,  1, M1, M2);
  AmpV.m_V[0][1] = VGamma(p1, p2, k, sigma,  1, -1, M1, M2);
  AmpV.m_V[1][0] = VGamma(p1, p2, k, sigma, -1,  1, M1, M2);
  AmpV.m_V[1][1] = VGamma(p1, p2, k, sigma, -1, -1, M1, M2);

}


Complex Ceex_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j,
                          double M1, double M2) {
  // eq 222 hep-ph/0006359
  return UGamma(p1, p2, k, h1, -i, -j, true, M1, M2);
}


Complex Ceex_Base::Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel) {
  Complex s, b1, b2;
  if (Sminus(k, p1) != -conj(Splus(k, p1))) {
    msg_Error() << "Wrong soft factors in " << METHOD << std::endl;
  }
  if (hel == -1) {
    b1 = sqrt(2) * Xi(p1, k) * Sminus(k, p1);
    b2 = sqrt(2) * Xi(p2, k) * Sminus(k, p2);
  }
  else {
    b1 = sqrt(2) * Xi(p1, k) * Splus(k, p1);
    b2 = sqrt(2) * Xi(p2, k) * Splus(k, p2);
  }
  m_b1 = Splus(k, p1);
  m_b2 = Splus(k, p2);
  m_Soft = -0.5 * (b1 / (p1 * k) - b2 / (p2 * k));
  s = -0.5 * (b1 / (p1 * k) - b2 / (p2 * k));
  m_pp1 = p1;
  m_pp2 = p2;
  m_kk = k;
  return -(m_e * s);
}
