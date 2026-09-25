/*!
  \file Ceex_Spinors.C

  The GPS spinor algebra of hep-ph/0006359: the T/U structures, the basic
  spinor products, the photon emission matrices U/V, and the soft factor.

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
#include "ATOOLS/Phys/Spinor.H"
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


/*!
  One leg's contribution to the eikonal current, b(p,k)/(p.k), normalised so
  that a charge-neutral pair reproduces Sfactor exactly:

      Sfactor(p1,p2,k,hel)  ==  SfactorLeg(p1,k,hel) - SfactorLeg(p2,k,hel)

  Sfactor computes that difference in one expression and is kept for the
  two-leg case; a stage with more than two legs needs the per-leg terms
  summed against their own weights w = Q*theta, which is what this exposes.
  The two routes differ only in the order the rounding happens.
*/
/*
  The sign of the negative-helicity branch.

  Sminus(k,p) = -conj(Splus(k,p)) - the check in Sfactor below asserts it - so
  writing s(-) with Sminus rather than conj(Splus) costs a relative minus sign
  between the two photon helicities. On its own that is a convention: flipping
  s(k,-1) everywhere multiplies each photon configuration by a global
  (-1)^#(sigma=-1), and |M|^2 is unchanged. It stops being free the moment the
  eikonal is subtracted from an amplitude computed elsewhere, because that
  amplitude carries ITS convention for epsilon_-.

  Comix is now the source of every Born and real amplitude here, and it uses
  the opposite sign. Measured, not assumed: with conserving kinematics the
  soft limit M_1 -> s beta_0 comes out at Re<M_1|s beta_0>/|M_1||s beta_0| =
  +1 for hel=+1 and -1 for hel=-1 (97 of 100 events). The factor of `hel`
  below is that measurement. It is invisible to any |s|^2 or sum-over-helicity
  test, which is why it survived this long.
*/
Complex Ceex_Base::SfactorLegGPS(const Vec4D &p, const Vec4D &k, int hel) {
  const Complex b(sqrt(2) * Xi(p, k)
                  * (hel == -1 ? -Sminus(k, p) : Splus(k, p)));
  return 0.5 * m_e * b / (p * k);
}

Vec4C Ceex_Base::ComixPolarisation(const Vec4D &k, int hel) const
{
  /*
    A transcription of METOOLS::CV::ConstructJ for a massless outgoing
    vector: EP for the current labelled H = 0, EM for H = 1, then the
    complex conjugate because the leg is outgoing. The gauge vector is the
    one Amplitude::SetGauge(0) leaves installed after the gauge test. The
    Spinor gauge permutation (R1,R2,R3) is whatever the Spinor class holds
    at the time, which is what Comix uses at the same time.
  */
  typedef ATOOLS::Spinor<double> Sp;
  static const Vec4D kg(1.0, 0.0, 1.0, 0.0);
  const Sp kp(1, kg), km(-1, kg);
  Vec4D p(k);
  if (p[1] == 0. && p[2] == 0.)   // the on-axis fix-up ConstructJ applies
    p[0] = p[0] < 0. ? -std::abs(p[3]) : std::abs(p[3]);
  const int nlg(Amplitude::s_nlegs);
  const int pflip(m_comixflip >= 0 ? (m_comixflip >> nlg) & 1
                                   : (m_comixphoflip & 1));
  const int hgi((hel > 0 ? 0 : 1) ^ pflip);
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
  One leg's contribution to the eikonal current, e (p.eps*)/(p.k), with the
  polarisation vector taken from Comix's own construction so that the
  subtraction M_1 - s B is a difference of like conventions. A stage's soft
  factor is the sum over its legs against w = Q*theta; a charge-neutral pair
  reproduces Sfactor exactly, which is now DEFINED as that difference.

  The GPS form this replaces (SfactorLegGPS) is the same vector up to a phase
  that depends on the photon's azimuth about the reference axes - trivially
  zero for a photon along the beams, which is why every check made on
  beam-collinear ISR photons passed, and 130 degrees on the wide-angle FSR
  photons of the n = 2 point that KKMC puts at real/Born = -0.17 and we had at
  +5.3. Common to every leg, so |s|^2, the I/F relative phase and hence rho_0
  are all unchanged by the switch.
*/
Complex Ceex_Base::SfactorLeg(const Vec4D &p, const Vec4D &k, int hel) {
  const Vec4C e(ComixPolarisation(k, hel));
  const Complex pe(p[0]*e[0] - p[1]*e[1] - p[2]*e[2] - p[3]*e[3]);
  return m_sfacphase * m_e * pe / (p * k);
}


Complex Ceex_Base::Sfactor(const Vec4D &p1, const Vec4D &p2, const Vec4D &k, int hel) {
  return SfactorLeg(p1, k, hel) - SfactorLeg(p2, k, hel);
}
