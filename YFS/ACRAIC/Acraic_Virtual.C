/*!
  \file Acraic_Virtual.C

  beta_0^1: the O(alpha) virtual. ISR/FSR vertex form factors (CnuA, exact
  in the masses) and the gamma-gamma / gamma-Z boxes (virtual IFI).
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "YFS/Main/Beta0_2.H"
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




Complex Acraic_Base::BoxGG() {
  const Complex MG(1e-30, 0.);
  const Complex L(log(-m_Tc / m_Sc));
  const Complex t1(log(m_Tc / m_Uc) * (log(MG * MG / m_Sc) + m_I * M_PI));
  const Complex t2(0.5 * m_Sc * (m_Uc - m_Tc) / sqr(m_Uc)
                   * (0.5 * L * L + m_I * M_PI * L));
  const Complex t3(-0.5 * m_Sc / m_Uc * (L + m_I * M_PI));
  return t1 + t2 + t3;
}


Complex Acraic_Base::BoxGZ() {
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



Complex Acraic_Base::BoxSubtract() {
  Complex MG = 1e-30;
  Complex sub = log(m_Tc / m_Uc) * log(MG * MG / sqrt(m_Uc * m_Tc))
                + 0.5 * log(m_Tc / m_Uc);
  return sub;
}


void Acraic_Base::MakeBoxMandelstams(const Vec4D &PX)
{
  if (m_pacraic.size() < 4) return;
  const double S(PX.Abs2());
  if (S <= 0.) return;
  Vec4D pd(m_pacraic[0] - m_pacraic[1]), qd(m_pacraic[2] - m_pacraic[3]);
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


Complex Acraic_Base::CnuA(double svar, double m1, double m2) const
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


void Acraic_Base::InfraredSubtractedME_0_1() {
  // m_Sc/m_Tc/m_Uc are set per partition by the caller from X = P - sum(ISR k).
  
  static const bool useboxessetting(Settings::GetMainSettings()["ACRAIC"]["BOXES"]
                                    .SetDefault(true).Get<bool>());
  // the boxes are the virtual IFI: none with YFS: IFI false or ISR only
  const bool useboxes(useboxessetting && m_ifi && !m_isronly);

  // The ut boxes are the tu ones with t and u swapped.
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
    // ACRAIC: RESONANCE_EXPONENTIATION: the resonant log sits in every Z
    // propagator, so remove it from the box (KKMC KeyINT 2, IntReson).
    Complex intreson(0., 0.);
    if (ResonanceExponentiation()) {
      const Complex mb2(sqr(m_MZ), -m_MZ * m_gZ);
      intreson = coef * (-2. * log(m_Tc / m_Uc) * log((mb2 - m_Sc) / mb2));
    }
    m_BoxGZtu -= intreson;
    const Complex t1(m_Tc), u1(m_Uc);
    m_Tc = u1; m_Uc = t1;
    m_BoxGGut = coef * (-BoxGG()) - SubBox;
    m_BoxGZut = coef * (-BoxGZ()) - SubBox - intreson;
    m_Tc = t1; m_Uc = u1;
  }

  
  const Complex LE(CnuA(m_s,  m_mass_I, m_mass_I) - 1.);
  const Complex LF(CnuA(m_sQ, m_mass_F, m_mass_F) - 1.);
  const Complex deltI(0.5 * sqr(m_qe) * m_alpha / Complex(M_PI, 0) * LE);
  const Complex deltF(HasFSR() ? 0.5 * sqr(m_qf) * m_alpha / Complex(M_PI, 0) * LF
                               : Complex(0., 0.));
  m_vertexI = deltI;
  m_vertexF = deltF;
  // Bhabha t-channel: both exchange vertices take the form factor at t
  const Complex deltT(m_bhabha
      ? 0.5 * sqr(m_qe) * m_alpha / Complex(M_PI, 0)
        * (CnuA(m_tinv, m_mass_I, m_mass_I) - 1.)
      : Complex(0., 0.));

  Amplitude born;
  BornAmplitude(m_pacraic, born, -1., -1., 0);

  Complex TG[2], TZ[2], UG[2], UZ[2];
  for (int j = 0; j <= 1; ++j) {
    const double h(1. - 2.*j);
    TG[j] = CouplingG() * m_propG;  TZ[j] = CouplingZ(h, 1) * m_propZ;
    UG[j] = CouplingG() * m_propG;  UZ[j] = CouplingZ(h, 0) * m_propZ;
  }

  /*
    ACRAIC: BETA0_2: two-loop vertex factors (KKMC MakeF1ini/MakeF1fin),
    F1_2 = F1_1 + F1_1^2/2 + 2 c_NLL (alpha/pi) Q^2 F1_1 (llnll), applied as
    (1 + F_ini)(1 + F_fin). The boxes are not dressed.
  */
  static const beta02::code b02(Settings::GetMainSettings()["ACRAIC"]["BETA0_2"]
                                .Get<beta02::code>());
  const Complex FI(b02 == beta02::none ? deltI
                   : Beta02VertexFactor(deltI, m_alpha, sqr(m_qe), b02));
  const Complex FF(b02 == beta02::none || !HasFSR() ? deltF
                   : Beta02VertexFactor(deltF, m_alpha, sqr(m_qf), b02));
  const Complex FT(b02 == beta02::none || !m_bhabha ? deltT
                   : Beta02VertexFactor(deltT, m_alpha, sqr(m_qe), b02));
  const Complex vert((1. + FI) * (1. + FF) - 1.);
  // Bhabha: s-channel Born dressed at s and sQ, t-channel Born at t
  const Complex vertT((1. + FT) * (1. + FT) - 1.);
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
          // Split the Born into its s- and t-channel parts (as in BornAmplitude).
          const Complex mt(m_bhabha
              ? m_TCt[j1]*m_Tampt[j1][j2][j3][j4]
              + m_UCt[j1]*m_Uampt[j1][j2][j3][j4]
              : Complex(0., 0.));
          const Complex ms(born.m_A[Idx(j1,j2,j3,j4)] - mt);
          const Complex corr(fac * (ms*vert + mt*vertT + boxy));
          m_AmpExpo1.m_A[Idx(j1,j2,j3,j4)] += corr;
          m_AmpBornVirt.m_A[Idx(j1,j2,j3,j4)] += corr;
          m_snapVirt.m_A[Idx(j1,j2,j3,j4)] = corr;   // the beta_0^1 increment
          m_beta01 += corr;   // scalar, for the diagnostics only
          // V = 1 + corr/born per partition and helicity, for the spread check
          // (both sides must be this partition's values)
          if (m_vpon) {
            const int f(Idx(j1,j2,j3,j4));
            if ((int)m_vpmin.size() <= f) { m_vpmin.resize(f+1, 1e30);
                                            m_vpmax.resize(f+1, -1e30); }
            // live helicities only: in mass-suppressed ones V is an artefact
            double bmax(0.);
            for (int g = 0; g < Amplitude::NHel(); ++g)
              bmax = Max(bmax, std::abs(born.m_A[g]));
            // divide by fac*born, since corr carries the same soft-factor fac
            const Complex b0(fac*born.m_A[f]);
            if (std::abs(b0) > 1e-3*std::abs(fac)*bmax) {
              const double v(std::abs(1. + corr/b0));
              m_vpmin[f] = Min(m_vpmin[f], v);
              m_vpmax[f] = Max(m_vpmax[f], v);
            }
          }
        }

}
