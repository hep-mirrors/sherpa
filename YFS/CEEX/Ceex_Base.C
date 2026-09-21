/*!
  \file Ceex_Base.C

  Construction, configuration and the per-event scaffolding: couplings,
  propagators, the CEEX momentum set, and MakeRho.

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

Amplitude::Amplitude() {
  for (int h0 = 0; h0 <= 1; ++h0) {
    for (int h1 = 0; h1 <= 1; ++h1) {
      for (int h2 = 0; h2 <= 1; ++h2) {
        for (int h3 = 0; h3 <= 1; ++h3) {
          m_A[h0][h1][h2][h3] = Complex(0, 0);
        }
      }
    }
  }
}


Ceex_Base::Ceex_Base(const Flavour_Vector &flavs)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  Settings& ss = Settings::GetMainSettings();
  m_onlyz = s["ONLYZ"].Get<int>();
  m_onlyg = s["ONLYG"].Get<int>();
  m_checkxs = s["CHECK_XS"].Get<int>();
  string widthscheme = ss["WIDTH_SCHEME"].Get<string>();
  // Both of Sherpa's width schemes put a CONSTANT M*Gamma in the propagator:
  // "Fixed" by definition, and "CMS" because the complex mass
  // mu^2 = M^2 - i*M*Gamma gives 1/(s - M^2 + i*M*Gamma). CEEX already takes a
  // complex sin^2(theta_W) from CMS, so the propagator has to follow it.
  //
  // KKMC instead uses an s-DEPENDENT width by default (CEEX.f:877-883; the
  // constant form is reachable there only via KeyZet = -1). Following KKMC here
  // disagreed with Sherpa's own Comix by up to +2.8% in |propZ|^2 at
  // x = Gamma_Z/M_Z = 0.027, i.e. right where radiative return crosses the Z
  // width. To compare against KKMC, set WIDTH_SCHEME: Fixed here AND run KKMC
  // with KeyZet = -1, so that both sides use the constant width.
  m_fixedwidth = (widthscheme == "Fixed" || widthscheme == "CMS");
  m_flavs = flavs;
  if (flavs.size() != 4) {
    THROW(fatal_error, "CEEX is only for 2->2");
  }
  if (flavs[2].IsNeutrino() && flavs[3].IsNeutrino()) {
    m_onlyz = true;
  }

  m_Q1Q2I = flavs[0].Charge() * flavs[1].Charge();
  m_QIQF  = flavs[0].Charge() * flavs[2].Charge();
  // Bhabha: the final pair is the beam pair, so a t-channel gamma/Z is exchanged
  // between the two fermion lines on top of the s-channel annihilation.
  m_bhabha = (flavs[0] == flavs[2]);
  m_Q1Q2F = flavs[2].Charge() * flavs[3].Charge();
  m_MZ = Flavour(kf_Z).Mass();
  m_gZ = Flavour(kf_Z).Width();
  double mw = Flavour(kf_Wplus).Mass();
  double MH = Flavour(kf_h0).Mass();
  double  GH  = Flavour(kf_h0).Width();
  double  GW  = Flavour(kf_Wplus).Width();
  double  GZ  = Flavour(kf_Z).Width();
  m_I   = Complex(0., 1.);

  double F_L = 0.;
  double F_R = 0.;

  m_sin2tw = MODEL::s_model->ComplexConstant("csin2_thetaW");
  if (Settings::GetMainSettings()["CEEX"]["REAL_SIN2THETAW"]
      .SetDefault(0).Get<int>() != 0)
    m_sin2tw = Complex(m_sin2tw.real(), 0.);
  m_sW = m_sin2tw;
  m_e = sqrt(4.*M_PI * m_alpha);
  m_cW = 1. - m_sW;
  m_norm = sqrt(16. * m_sW * (1. - m_sW));
  m_qe       = m_flavs[0].Charge();
  m_qf       = m_flavs[2].Charge();
  m_Q1Q2I = m_flavs[0].Charge() * m_flavs[1].Charge();
  m_ae       = 2.*m_flavs[0].IsoWeak();
  m_af       = 2.*m_flavs[2].IsoWeak();
  // Keep 2*T3 and 4*Q*sw^2 separately: the electroweak kappa factors multiply
  // only the sin^2 piece (GPS_EWFFact), so the two cannot be pre-combined.
  m_t3e2 = m_ae;
  m_t3f2 = m_af;
  m_qesw = 4.*m_qe * m_sin2tw;
  m_qfsw = 4.*m_qf * m_sin2tw;
  m_ve       = (m_ae - 4.*m_qe * m_sin2tw) / m_norm;
  m_vf       = (m_af - 4.*m_qf * m_sin2tw) / m_norm;
  m_ae /= m_norm;
  m_af /= m_norm;
  m_weak = s["WEAK"].Get<int>();
  m_mass_I = flavs[0].Mass();
  m_mass_F = flavs[2].Mass();
  // full EW couplings
  m_I_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[0].IsoWeak()
          - 2.*flavs[0].Charge() * m_sW * m_sW);

  m_I_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[0].Charge() * m_sW * m_sW);

  m_F_L = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*flavs[2].IsoWeak()
          - 2.*flavs[2].Charge() * m_sW * m_sW);

  m_F_R = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*flavs[2].Charge() * m_sW * m_sW);
  m_cL = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (2.*m_flavs[1].IsoWeak()
         - 2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_cR = -m_I * sqrt(4 * M_PI * m_alpha) / (2.*m_sW * m_cW) * (-2.*m_flavs[1].Charge() * m_sW * m_sW) / m_norm;
  m_zeta = {1, 1, 0, 0};
  m_eta = {0, 0, 1, 0};
  m_b  = {0.0,  0.8723e0, -0.7683e0, 0.3348e0};
}


void Ceex_Base::RegisterDefaults()
{
  Scoped_Settings s{ Settings::GetMainSettings()["CEEX"] };
  s["ONLYZ"].SetDefault(0);
  s["ONLYG"].SetDefault(0);
  s["CHECK_XS"].SetDefault(0);
  s["WEAK"].SetDefault(1);
}



void Ceex_Base::Init(const Vec4D_Vector &p)
{
  m_momenta = p;
  m_cms = Poincare(m_momenta[0] + m_momenta[1]);
  Poincare cms = m_cms;
  Poincare Rot = Poincare(Vec4D(0., 0., 0., 1.));
  for (size_t i(0); i < p.size(); ++i) {
    cms.Boost(m_momenta[i]);
    // cms.Boost(m_bornmomenta[i]);
  }
  m_crude = 2.0 / (4.0 * M_PI);
  for (size_t i(0); i < m_isrphotons.size(); ++i) {
    m_crude /= pow(2 * M_PI, 3);
  }
  if (m_momenta.size() >= 6) {
    m_sp = (m_momenta[4] + m_momenta[5]).Abs2();
    m_sQ = m_sp;
  } else {
    m_sp = (m_momenta[2] + m_momenta[3]).Abs2();
    m_sQ = m_sp;
  }
  m_T = 0;
}



void Ceex_Base::MakeProp()
{
  if (m_fixedwidth) {
    m_propZ =  1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_MZ);
  }
  else {
    m_propZ =   1. / Complex(m_sp - sqr(m_MZ), m_gZ * m_sp / m_MZ);
  }
  m_propG = 1. / Complex(m_sp,0);
  if (m_onlyz)  m_propG = 0;
  if (m_onlyg)  m_propZ = 0;
  m_prop =  m_propZ + m_propG;
}


void Ceex_Base::MakePropT(const Vec4D_Vector &p)
{
  if (!m_bhabha) {
    m_propGt = m_propZt = Complex(0., 0.);
    return;
  }
  m_tinv = (p[0] - p[2]).Abs2();
  // Fixed width, not the running form: t is spacelike, so a width scaled by t
  // would put the pole on the wrong side. The width is kept rather than dropped
  // because the complex-mass scheme carries the complex mass in every
  // propagator, spacelike included - that is what keeps it gauge invariant, and
  // it is the same scheme the complex sin^2(theta_W) above comes from. Dropping
  // it costs a factor 12 in the agreement with Comix across scattering angle.
  m_propZt = 1. / Complex(m_tinv - sqr(m_MZ), m_gZ * m_MZ);
  m_propGt = 1. / Complex(m_tinv, 0.);
  if (m_onlyz) m_propGt = 0.;
  if (m_onlyg) m_propZt = 0.;
}


void Ceex_Base::MakeEWFF(double svar, double costhd)
{
  m_kapE = m_kapF = m_kapEF = Complex(1., 0.);
  m_rhoEW = m_gamVPi = Complex(1., 0.);
  m_vvcor = Complex(1., 0.);
  if (!m_weak) {           // tree couplings == KKMC with KeyElw = 0
    m_ve = (m_t3e2 - m_qesw) / m_norm;
    m_vf = (m_t3f2 - m_qfsw) / m_norm;
    return;
  }
  // --- weak form factors go here; not yet implemented ---
  m_ve = (m_t3e2 - m_qesw * m_kapE) / m_norm;
  m_vf = (m_t3f2 - m_qfsw * m_kapF) / m_norm;
  // Angle-dependent double-vector correction; kapEF carries the box content.
  const Complex vvcef((m_t3e2*m_t3f2
                       - m_qesw*m_t3f2*m_kapE
                       - m_qfsw*m_t3e2*m_kapF
                       + m_qesw*m_qfsw*m_kapEF) / (m_norm*m_norm));
  m_vvcor = (std::abs(m_ve*m_vf) > 0.) ? vvcef/(m_ve*m_vf) : Complex(1., 0.);
}


Complex Ceex_Base::CouplingZ(double  j, int mode) {
  if (m_onlyg) return 0.;
  Complex zcpl;
  if (mode == 1) {
    zcpl = m_ve * m_vf * m_vvcor - dcmplx(j) * m_ae * m_vf + dcmplx(j) * m_ve * m_af - m_ae * m_af;
  }
  else if (mode == 0) {
    zcpl = m_ve * m_vf * m_vvcor - dcmplx(j) * m_ae * m_vf - dcmplx(j) * m_ve * m_af + m_af * m_ae;
  }
  else msg_Error() << METHOD << "\n wrong mode\n";

  if (zcpl == 0.) {
    msg_Error() << "Z coupling is Zero!\n";
  }
  return zcpl;
}



Complex Ceex_Base::CouplingG() {
  m_gcpl = Complex(m_QIQF, 0);
  return m_gcpl;
}


void Ceex_Base::BuildCeexMomenta()
{
  m_pceex.clear();
  m_pceex.push_back(m_momenta[0]);
  m_pceex.push_back(m_momenta[1]);
  if (m_momenta.size() >= 6) {
    m_pceex.push_back(m_momenta[4]);   // physical outgoing fermion
    m_pceex.push_back(m_momenta[5]);
  } else {
    m_pceex.push_back(m_momenta[2]);
    m_pceex.push_back(m_momenta[3]);
  }
}


void Ceex_Base::ZerAmplit() {
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          m_AmpExpo0.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_AmpExpo1.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_AmpBornVirt.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_AmpBornReal.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_snapBorn.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_snapVirt.m_A[j1][j2][j3][j4] = Complex(0., 0.);
          m_snapReal.m_A[j1][j2][j3][j4] = Complex(0., 0.);
        }
}



void Ceex_Base::MakeRho() {
  double sum0(0.), sum1(0.);
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          sum0 += std::real(m_AmpExpo0.m_A[j1][j2][j3][j4]
                            * conj(m_AmpExpo0.m_A[j1][j2][j3][j4]));
          sum1 += std::real(m_AmpExpo1.m_A[j1][j2][j3][j4]
                            * conj(m_AmpExpo1.m_A[j1][j2][j3][j4]));
        }
  // Average over the four initial-state helicity configurations.
  m_result0 = sum0 / 4.;
  m_result  = sum1 / 4.;
  double sumbv(0.), sumbr(0.);
  for (int j1 = 0; j1 <= 1; ++j1)
    for (int j2 = 0; j2 <= 1; ++j2)
      for (int j3 = 0; j3 <= 1; ++j3)
        for (int j4 = 0; j4 <= 1; ++j4) {
          sumbv += std::real(m_AmpBornVirt.m_A[j1][j2][j3][j4]
                             * conj(m_AmpBornVirt.m_A[j1][j2][j3][j4]));
          sumbr += std::real(m_AmpBornReal.m_A[j1][j2][j3][j4]
                             * conj(m_AmpBornReal.m_A[j1][j2][j3][j4]));
        }
  m_resultbv = sumbv / 4.;
  m_resultbr = sumbr / 4.;
  m_rho0sum += m_result0;
  m_rho1sum += m_result;
  m_rhobvsum += m_resultbv;
  m_rhobrsum += m_resultbr;
  m_rhocrudsum += m_rhocrud;
  ++m_rhon;
}


void Ceex_Base::Reset() {
  m_result = 0;
}


double Ceex_Base::Xi(const Vec4D p, const Vec4D q) {
  return sqrt((m_zeta * p) / (q * m_zeta));
}

double Ceex_Base::RealFactorPhoton(size_t j) const
{
  if (j >= m_realphot.size() || m_result0 <= 0.) return 0.;
  double sum(0.);
  for (int a = 0; a <= 1; ++a)
    for (int b = 0; b <= 1; ++b)
      for (int c = 0; c <= 1; ++c)
        for (int d = 0; d <= 1; ++d) {
          const Complex z(m_AmpExpo0.m_A[a][b][c][d] + m_realphot[j].m_A[a][b][c][d]);
          sum += std::real(z * conj(z));
        }
  return sum/4./m_result0 - 1.;
}
