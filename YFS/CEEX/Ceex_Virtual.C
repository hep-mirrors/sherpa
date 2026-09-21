/*!
  \file Ceex_Virtual.C

  beta_0^1: the O(alpha) virtual. The initial- and final-state vertex form
  factors (CnuA, exact in the masses) and the gamma-gamma / gamma-Z boxes
  that carry the virtual initial-final interference.

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




Complex Ceex_Base::BoxGG() {
  const Complex MG(1e-30, 0.);
  const Complex L(log(-m_Tc / m_Sc));
  const Complex t1(log(m_Tc / m_Uc) * (log(MG * MG / m_Sc) + m_I * M_PI));
  const Complex t2(0.5 * m_Sc * (m_Uc - m_Tc) / sqr(m_Uc)
                   * (0.5 * L * L + m_I * M_PI * L));
  const Complex t3(-0.5 * m_Sc / m_Uc * (L + m_I * M_PI));
  return t1 + t2 + t3;
}


Complex Ceex_Base::BoxGZ() {
  Complex MG = 1e-30;
  Complex mb2 = Complex(sqr(m_MZ), -m_MZ * m_gZ);
  Complex t1 = log(m_Tc / m_Uc) * (log(MG * MG / sqrt(m_Tc * m_Uc)));
  Complex t2 = -2.*log(m_Tc / m_Uc) * log((mb2 - m_Sc) / mb2)
               + DiLog((mb2 + m_Uc) / mb2) - DiLog((mb2 + m_Tc) / mb2);

  Complex t3 = (mb2 - m_Sc) * (m_Uc - m_Tc - mb2) / (sqr(m_Uc)) * (
                 log(-m_Tc / m_Sc) * log((mb2 - m_Sc) / mb2)
                 + DiLog((mb2 + m_Tc) / mb2) - DiLog((mb2 - m_Sc) / mb2));
  Complex t4 = sqr(mb2 - m_Sc) / (m_Sc * m_Uc) * log((mb2 - m_Sc) / mb2)
               + (mb2 - m_Sc) / m_Uc * log(-m_Tc / mb2);

  return t1 + t2 + t3 + t4;

}



Complex Ceex_Base::BoxSubtract() {
  Complex MG = 1e-30;
  Complex sub = log(m_Tc / m_Uc) * log(MG * MG / sqrt(m_Uc * m_Tc))
                + 0.5 * log(m_Tc / m_Uc);
  return sub;
}


void Ceex_Base::MakeBoxMandelstams(const Vec4D &PX)
{
  if (m_pceex.size() < 4) return;
  const double S(PX.Abs2());
  if (S <= 0.) return;
  Vec4D pd(m_pceex[0] - m_pceex[1]), qd(m_pceex[2] - m_pceex[3]);
  pd = pd - ((PX*pd)/S) * PX;
  qd = qd - ((PX*qd)/S) * PX;
  const double a(pd.Abs2()), b(qd.Abs2());
  // KKMC's ThetaD treats a*b <= 0 as an error.
  if (a*b <= 0.) return;
  double cth(-(qd*pd)/sqrt(std::abs(a*b)));
  if (cth >  1.) cth =  1.;
  if (cth < -1.) cth = -1.;
  m_Sc = Complex(S, 0.);
  m_Tc = Complex(-S*(1. - cth)/2., 0.);
  m_Uc = Complex(-S*(1. + cth)/2., 0.);
}


Complex Ceex_Base::CnuA(double svar, double m1, double m2) const
{
  const double mm(m1*m2);
  if (std::abs(mm) < 1e-10)
    THROW(fatal_error, "CnuA needs two non-zero masses; "
                       "set Massive: true on both legs.");
  const Complex Mas12(mm, 0.), Nu((-svar + m1*m1 + m2*m2)/2., 0.);
  const Complex xlam(sqrt((Nu - Mas12)*(Nu + Mas12)));
  const Complex z(svar > 0. ? Mas12/(Nu - xlam) : (Nu + xlam)/Mas12);
  // KKbvir::CDLN(z, +1): log(z) except on the negative real axis, where the
  // +i eps prescription adds +i pi.
  const Complex lz((z.imag() == 0. && z.real() <= 0.)
                   ? log(-z) + Complex(0., M_PI) : log(z));
  return Nu/xlam * lz;
}


void Ceex_Base::InfraredSubtractedME_0_1() {
  // m_Sc/m_Tc/m_Uc are set per PARTITION by the caller, from that partition's
  // X = P - sum(ISR k) - they are not event constants.
  
  static const bool useboxes(Settings::GetMainSettings()["CEEX"]["BOXES"]
                             .SetDefault(1).Get<int>() != 0);

  // Boxes first. The ut combination is the tu one with t and u swapped.
  const Complex coef(m_alpha * m_qe * m_qf / M_PI);
  
  
  m_BoxGGtu = m_BoxGZtu = m_BoxGGut = m_BoxGZut = Complex(0., 0.);
  if (useboxes || m_checkxs) {
    const Complex SubBox(coef * BoxSubtract());
    if (m_checkxs) {
      m_rawboxGG = BoxGG(); m_rawboxGZ = BoxGZ(); m_rawboxSub = BoxSubtract();
      m_rawcoef  = coef;
    }
    m_BoxGGtu = coef * BoxGG() - SubBox;
    m_BoxGZtu = coef * BoxGZ() - SubBox;
    const Complex t1(m_Tc), u1(m_Uc);
    m_Tc = u1; m_Uc = t1;
    m_BoxGGut = coef * (-BoxGG()) - SubBox;
    m_BoxGZut = coef * (-BoxGZ()) - SubBox;
    m_Tc = t1; m_Uc = u1;
  }

  
  const Complex LE(CnuA(m_s,  m_mass_I, m_mass_I) - 1.);
  const Complex LF(CnuA(m_sQ, m_mass_F, m_mass_F) - 1.);
  const Complex deltI(0.5 * sqr(m_qe) * m_alpha / Complex(M_PI, 0) * LE);
  const Complex deltF(HasFSR() ? 0.5 * sqr(m_qf) * m_alpha / Complex(M_PI, 0) * LF
                               : Complex(0., 0.));
  m_vertexI = deltI;
  m_vertexF = deltF;
  // Bhabha t-channel vertex. Both vertices of the exchange diagram sit on an
  // electron line at spacelike momentum transfer t - the (p1,p3) line and the
  // (p2,p4) line - so each takes the SAME form factor evaluated at t. CnuA is
  // valid there: KKMC's BVR_CnuA header states it is "appropriate for s and
  // t-chanels" with no small-mass approximation, and the port keeps its
  // sign-of-svar stability branch.
  const Complex deltT(m_bhabha
      ? 0.5 * sqr(m_qe) * m_alpha / Complex(M_PI, 0)
        * (CnuA(m_tinv, m_mass_I, m_mass_I) - 1.)
      : Complex(0., 0.));

  Amplitude born;
  BornAmplitude(m_pceex, born, -1., -1., 0);

  Complex TG[2], TZ[2], UG[2], UZ[2];
  for (int j = 0; j <= 1; ++j) {
    const double h(1. - 2.*j);
    TG[j] = CouplingG() * m_propG;  TZ[j] = CouplingZ(h, 1) * m_propZ;
    UG[j] = CouplingG() * m_propG;  UZ[j] = CouplingZ(h, 0) * m_propZ;
  }

  const Complex vert((1. + deltI) * (1. + deltF) - 1.);
  // The s-channel factorisation does not carry over to Bhabha: the s-channel
  // Born takes its corrections at s and sQ, the t-channel Born takes both at t.
  // So the two Born pieces must be dressed separately rather than one
  // multiplicative factor being applied to their sum.
  const Complex vertT((1. + deltT) * (1. + deltT) - 1.);
  const Complex fac(m_e * m_e * m_cfac);
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          const int H1(1 - 2*j1), H3(1 - 2*j3);
          const Complex bgg((H1*H3 == 1) ? m_BoxGGtu : m_BoxGGut);
          const Complex bgz((H1*H3 == 1) ? m_BoxGZtu : m_BoxGZut);
          const Complex boxy(useboxes
                           ? m_Tamp[j1][j2][j3][j4]*(TG[j1]*bgg + TZ[j1]*bgz)
                           + m_Uamp[j1][j2][j3][j4]*(UG[j1]*bgg + UZ[j1]*bgz)
                           : Complex(0., 0.));
          // Split the Born back into its s- and t-channel parts (the same
          // decomposition BornAmplitude assembled) so each gets its own vertex.
          const Complex mt(m_bhabha
              ? m_TCt[j1]*m_Tampt[j1][j2][j3][j4]
              + m_UCt[j1]*m_Uampt[j1][j2][j3][j4]
              : Complex(0., 0.));
          const Complex ms(born.m_A[j1][j2][j3][j4] - mt);
          const Complex corr(fac * (ms*vert + mt*vertT + boxy));
          m_AmpExpo1.m_A[j1][j2][j3][j4] += corr;
          m_AmpBornVirt.m_A[j1][j2][j3][j4] += corr;
          m_snapVirt.m_A[j1][j2][j3][j4] = corr;   // the beta_0^1 increment
          m_beta01 += corr;   // scalar, for the diagnostics only
        }
}
