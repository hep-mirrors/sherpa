/*!
  \file Acraic_Spinors.C

  The GPS spinor algebra of hep-ph/0006359: the T/U structures, the basic
  spinor products, the photon emission matrices U/V, and the soft factor.

*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;



Complex Acraic_Base::T(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                     double M1, double M2) {
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = Xi(p1, p2);
  double sq2 = Xi(p2, p1);
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



Complex Acraic_Base::Tp(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
                      double M1, double M2) {
  h1 = -h1;
  h2 = -h2;
  Complex s(0, 0);
  const double m1(M1 >= 0. ? M1 : p1.Mass());
  const double m2(M2 >= 0. ? M2 : p2.Mass());
  double sq1 = Xi(p1, p2);
  double sq2 = Xi(p2, p1);
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


Complex Acraic_Base::U(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
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


Complex Acraic_Base::Up(const Vec4D &p1, const Vec4D &p2, int h1, int h2, int C1, int C2,
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


Complex Acraic_Base::S(const Vec4D &p1, const Vec4D &p2, int h1, int h2) {
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


Complex Acraic_Base::S(const Vec4D &p1, const Vec4D &p2, double m1, double m2, int h1, int h2) {
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



Complex Acraic_Base::Splus(const Vec4D &p, const Vec4D &q) {
  Complex sp = -Complex(q[2], q[3]) * sqrt((p[0] - p[1]) / (q[0] - q[1]));
  sp += Complex(p[2], p[3]) * sqrt((q[0] - q[1]) / (p[0] - p[1]));
  return sp;
}


Complex Acraic_Base::Sminus(const Vec4D &p, const Vec4D &q) {
  return -conj(Splus(p, q));
}



void Acraic_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpU,
                       double M1, double M2) {

  AmpU.m_U[0][0] = UGamma(p1, p2, k, sigma,  1,  1, false, M1, M2);
  AmpU.m_U[0][1] = UGamma(p1, p2, k, sigma,  1, -1, false, M1, M2);
  AmpU.m_U[1][0] = UGamma(p1, p2, k, sigma, -1,  1, false, M1, M2);
  AmpU.m_U[1][1] = UGamma(p1, p2, k, sigma, -1, -1, false, M1, M2);

}



Complex Acraic_Base::UGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j,
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
  Complex coeff1 = sqrt(2) * Xi(p1, k);
  Complex coeff2 = sqrt(2) * Xi(p2, k);
  if ( i == 1 && j == 1) return sqr2 * Xi(p2, k) * Splus(k, p1hat);
  else if (i == 1 && j == -1) return 0;
  else if (i == -1 && j == 1) return sqr2 * (m2 * Xi(p1, p2) - m1 * Xi(p2, p1));
  else if (i == -1 && j == -1) return sqr2 * Xi(p1, k) * Splus(k, p2hat);
  else {
    msg_Out() << "h1 = " << i << "\n h2 = " << j;
    THROW(fatal_error, "Wrong helicities in ACRAIC");
  }
}


void Acraic_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int sigma, Amplitude &AmpV,
                       double M1, double M2) {

  AmpV.m_V[0][0] = VGamma(p1, p2, k, sigma,  1,  1, M1, M2);
  AmpV.m_V[0][1] = VGamma(p1, p2, k, sigma,  1, -1, M1, M2);
  AmpV.m_V[1][0] = VGamma(p1, p2, k, sigma, -1,  1, M1, M2);
  AmpV.m_V[1][1] = VGamma(p1, p2, k, sigma, -1, -1, M1, M2);

}


Complex Acraic_Base::VGamma(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int h1, int i, int j,
                          double M1, double M2) {
  // eq 222 hep-ph/0006359
  return UGamma(p1, p2, k, h1, -i, -j, true, M1, M2);
}


//! One leg's eikonal term b(p,k)/(p.k) in the GPS photon convention.
// hel=-1 uses -Sminus = conj(Splus) to match Comix's epsilon_- sign.
Complex Acraic_Base::SfactorLegGPS(const Vec4D &p, const Vec4D &k, int hel) {
  const Complex b(sqrt(2) * Xi(p, k)
                  * (hel == -1 ? -Sminus(k, p) : Splus(k, p)));
  return 0.5 * m_e * b / (p * k);
}

Vec4C Acraic_Base::ComixPolarisation(const Vec4D &k, int hel) const
{
  // METOOLS::CV::ConstructJ for a massless outgoing vector (conjugated), with
  // the gauge vector Amplitude::SetGauge(0) leaves installed.
  typedef ATOOLS::Spinor<double> Sp;
  static const Vec4D kg(1.0, 0.0, 1.0, 0.0);
  const Sp kp(1, kg), km(-1, kg);
  Vec4D p(k);
  if (p[1] == 0. && p[2] == 0.)   // the on-axis fix-up ConstructJ applies
    p[0] = p[0] < 0. ? -std::abs(p[3]) : std::abs(p[3]);
  const int nlg(Amplitude::s_nlegs);
  const int pflip(m_comixflip >= 0 ? ComixPhotonFlip(nlg)
                                   : LegBitOf(m_comixphoflip, 0));
  const int hgi(FlipLegs(HelicityBit<int>(hel), pflip));
  Vec4C e;
  Complex den;
  if (hgi == 0) {                  // EP
    const Sp a(kp), b(-1, p);
    e[0] = a.U1()*b.U1() + a.U2()*b.U2();
    e[Sp::R3()] = a.U1()*b.U1() - a.U2()*b.U2();
    e[Sp::R1()] = a.U1()*b.U2() + a.U2()*b.U1();
    e[Sp::R2()] = Complex(0., 1.)*(a.U1()*b.U2() - a.U2()*b.U1());
    den = sqrt(2.) * std::conj(km*b);
  } else {                         // EM
    const Sp a(1, p), b(km);
    e[0] = a.U1()*b.U1() + a.U2()*b.U2();
    e[Sp::R3()] = a.U1()*b.U1() - a.U2()*b.U2();
    e[Sp::R1()] = a.U1()*b.U2() + a.U2()*b.U1();
    e[Sp::R2()] = Complex(0., 1.)*(a.U1()*b.U2() - a.U2()*b.U1());
    den = sqrt(2.) * std::conj(kp*a);
  }
  for (int i(0); i < 4; ++i) e[i] = std::conj(e[i]/den);
  return e;
}

/*!
  One leg's eikonal term e (p.eps*)/(p.k) with Comix's polarisation, so
  M_1 - s B subtracts like conventions; a stage sums its legs with w = Q*theta.
*/
Complex Acraic_Base::SfactorLeg(const Vec4D &p, const Vec4D &k, int hel) {
  const Vec4C e(ComixPolarisation(k, hel));
  const Complex pe(p[0]*e[0] - p[1]*e[1] - p[2]*e[2] - p[3]*e[3]);
  return m_sfacphase * m_e * pe / (p * k);
}

Complex Acraic_Base::SfactorLeg(const Vec4D &pnum, const Vec4D &ppole,
                              const Vec4D &k, int hel) {
  const Vec4C e(ComixPolarisation(k, hel));
  const Complex pe(pnum[0]*e[0] - pnum[1]*e[1] - pnum[2]*e[2] - pnum[3]*e[3]);
  return m_sfacphase * m_e * pe / (ppole * k);
}


Complex Acraic_Base::Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel) {
  return SfactorLeg(p1, k, hel) - SfactorLeg(p2, k, hel);
}
