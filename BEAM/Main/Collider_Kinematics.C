#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "BEAM/Main/Collider_Kinematics.H"

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Collider_Kinematics::Collider_Kinematics(Beam_Base **beams)
    : Kinematics_Base(beams), m_mode(collidermode::unknown) {
  if (p_beams[0]->Type() == beamspectrum::monochromatic &&
      p_beams[1]->Type() == beamspectrum::monochromatic)
    m_mode = collidermode::monochromatic;
  else if (p_beams[0]->Type() != beamspectrum::monochromatic &&
           p_beams[1]->Type() == beamspectrum::monochromatic)
    m_mode = collidermode::spectral_1;
  else if (p_beams[0]->Type() == beamspectrum::monochromatic &&
           p_beams[1]->Type() != beamspectrum::monochromatic)
    m_mode = collidermode::spectral_2;
  else if (p_beams[0]->Type() != beamspectrum::monochromatic &&
           p_beams[1]->Type() != beamspectrum::monochromatic)
    m_mode = collidermode::both_spectral;
  if (m_mode == collidermode::unknown)
    THROW(fatal_error, "Bad settings for collider mode.");
  InitSystem();
  InitIntegration();
}

Collider_Kinematics::~Collider_Kinematics() {}

void Collider_Kinematics::InitSystem() {
  m_Ecms = sqrt(m_S);

  rpa->gen.SetEcms(m_Ecms);
  rpa->gen.SetPBeam(0, p_beams[0]->InMomentum());
  rpa->gen.SetPBeam(1, p_beams[1]->InMomentum());

  Settings::GetMainSettings().AddGlobalTag("E_CMS", ToString(m_Ecms));
  m_on = (m_mode != collidermode::monochromatic);
  m_x[0] = 1.;
  m_x[1] = 1.;
  m_LabBoost = Poincare(p_beams[0]->InMomentum() + p_beams[1]->InMomentum());
  m_CMSBoost = Poincare(p_beams[0]->OutMomentum() + p_beams[1]->OutMomentum());
}

void Collider_Kinematics::InitIntegration() {
  Beam_Parameters parameters;
  // check for if they have been initialised to other values
  double sminratio = parameters("BEAM_SMIN");
  double smaxratio = parameters("BEAM_SMAX");
  m_xmin = p_beams[0]->Xmin() * p_beams[1]->Xmin();
  m_xmax = p_beams[0]->Xmax() * p_beams[1]->Xmax();
  m_smin = m_S * Max(m_xmin, sminratio);
  m_smax = m_S * Min(m_xmax, smaxratio);
  // TODO: the rapidity interval can be done in a better way.
  m_ymin = -20.;
  m_ymax = 20.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * (p_beams[0]->Exponent() + p_beams[1]->Exponent());
}

bool Collider_Kinematics::operator()(ATOOLS::Vec4D *moms) {
  m_sprime = m_sprimekey[3];
  if (m_sprime < m_sprimekey[0] || m_sprime > m_sprimekey[1]) {
    msg_Error() << METHOD << "(..): " << om::red << "s' out of bounds.\n"
                << om::reset
                << "  s'_{min}, s'_{max 1,2} vs. s': " << m_sprimekey[0] << ", "
                << m_sprimekey[1] << ", " << m_sprimekey[2] << " vs. "
                << m_sprime << std::endl;
    return false;
  }
  CalculateSudakovMomenta();
  switch (m_mode) {
  case collidermode::monochromatic:
    return MakeMonochromaticBeams();
  case collidermode::spectral_1:
    return MakeSpectral1Beams(moms);
  case collidermode::spectral_2:
    return MakeSpectral2Beams(moms);
  case collidermode::both_spectral:
    return MakeCollinearBeams(moms);
  case collidermode::unknown:
    THROW(fatal_error,
          "Unknown collider mode, impossible to build kinematics.");
  }
  return false;
}

bool Collider_Kinematics::MakeMonochromaticBeams() {
  // Should actually not be called, because the beamhandler is considered to be
  // off for monochromatic beams
  return true;
}

bool Collider_Kinematics::MakeCollinearBeams(ATOOLS::Vec4D *moms) {
  double tau = CalculateTau();
  double yt = exp(m_ykey[2] - 0.5 * log((tau + m_m2[1]) / (tau + m_m2[0])) -
                  m_Plab.Y());
  m_x[0] = m_xkey[4] = sqrt(tau) * yt;
  m_x[1] = m_xkey[5] = sqrt(tau) / yt;
  moms[0] = m_x[0] * m_p_plus + m_m2[0] / m_S / m_x[0] * m_p_minus;
  moms[1] = m_x[1] * m_p_minus + m_m2[1] / m_S / m_x[1] * m_p_plus;
  for (size_t i = 0; i < 2; ++i) {
    p_beams[i]->SetOutMomentum(moms[i]);
  }
  m_CMSBoost = Poincare(moms[0] + moms[1]);
  return true;
}

bool Collider_Kinematics::MakeSpectral1Beams(ATOOLS::Vec4D *moms) {
  m_x[1] = m_xkey[5] =
      1.; // Should actually be p_beams[1]->InMomentum().PMinus() /
          // m_p_minus.PMinus(), but leads to violation of momentum conservation
          // during event generation
  m_x[0] = m_xkey[4] = CalculateTau() / m_x[1];
  moms[0] = m_x[0] * m_p_plus + m_m2[0] / m_S / m_x[0] * m_p_minus;
  moms[1] = p_beams[1]->InMomentum();
  for (size_t i = 0; i < 2; ++i) {
    p_beams[i]->SetOutMomentum(moms[i]);
  }
  m_CMSBoost = Poincare(moms[0] + moms[1]);
  return true;
}

bool Collider_Kinematics::MakeSpectral2Beams(ATOOLS::Vec4D *moms) {
  m_x[0] = m_xkey[4] =
      1.; // Should actually be p_beams[0]->InMomentum().PPlus() /
          // m_p_plus.PPlus(), see above
  m_x[1] = m_xkey[5] = CalculateTau() / m_x[0];
  moms[0] = p_beams[0]->InMomentum();
  moms[1] = m_x[1] * m_p_minus + m_m2[1] / m_S / m_x[1] * m_p_plus;
  for (size_t i = 0; i < 2; ++i) {
    p_beams[i]->SetOutMomentum(moms[i]);
  }
  m_CMSBoost = Poincare(moms[0] + moms[1]);
  return true;
}

double Collider_Kinematics::CalculateTau() {
  double tau = (m_sprime - m_m2[0] - m_m2[1]) / m_S / 2.;
  if (tau * tau < m_m2[0] * m_m2[1] / (m_S * m_S)) {
    msg_Error() << METHOD << "(): s' out of range." << std::endl;
    return false;
  }
  tau += sqrt(tau * tau - m_m2[0] * m_m2[1] / (m_S * m_S));
  return tau;
}

void Collider_Kinematics::CalculateSudakovMomenta() {
  Vec4D pa = p_beams[0]->InMomentum();
  Vec4D pb = p_beams[1]->InMomentum();
  double gam = pa * pb + sqrt(sqr(pa * pb) - pa.Abs2() * pb.Abs2());
  double bet = 1.0 / (1.0 - pa.Abs2() / gam * pb.Abs2() / gam);
  m_p_plus = bet * (pa - pa.Abs2() / gam * pb);
  m_p_minus = bet * (pb - pb.Abs2() / gam * pa);
}

void Collider_Kinematics::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid + string("s'"), 5, 0, info);
  m_ykey.Assign(m_keyid + string("y"), 3, 0, info);
  // Convention for m_xkey:
  // [x_{min,beam0}, x_{min,beam1}, x_{max,beam0}, x_{max,beam1}, x_{val,beam0},
  // x_{val,beam1}] The limits, i.e. index 0,1,2,3 are saved as log(x), the
  // values are saved linearly.
  m_xkey.Assign(m_keyid + string("x"), 6, 0, info);
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_smax;
  m_sprimekey[2] = m_S;
  m_sprimekey[3] = m_S;
  m_sprimekey[4] = -m_S;
  m_ykey[0] = m_ymin;
  m_ykey[1] = m_ymax;
  m_ykey[2] = 0.;
}

void Collider_Kinematics::SetLimits() {
  m_sprimekey[0] = Max(m_smin, m_sminPS);
  m_sprimekey[1] = m_sprimekey[2] = m_smax;
  m_sprimekey[3] = m_S;
  m_ykey[0] = m_ymin;
  m_ykey[1] = m_ymax;
  m_ykey[2] = 0.;
  for (size_t i = 0; i < 2; i++) {
    double p = i == 0 ? p_beams[0]->OutMomentum().PPlus()
                      : p_beams[1]->OutMomentum().PMinus();
    double e = p_beams[i]->OutMomentum()[0];
    m_xkey[i] =
        (IsZero(m_m[i], 1.e-13) ? -0.5 * std::numeric_limits<double>::max()
                                : 2. * log(m_m[i] / p));
    m_xkey[i + 2] = log(
        Min(p_beams[i]->Xmax(), (e / p * (1.0 + sqrt(1.0 - sqr(m_m[i] / e))))));
    m_xkey[i + 4] = m_x[i];
  }
  // sprime's with masses - still need to check for masses
  double sprimemin = Max(m_sprimekey[0], m_S * exp(m_xkey[0] + m_xkey[1]));
  if (sprimemin > sqr(m_m[0] + m_m[1]))
    m_sprimekey[0] = sprimemin;
  double sprimemax = Min(m_smax, m_S * exp(m_xkey[2] + m_xkey[3]));
  if (sprimemax > sqr(m_m[0] + m_m[1]))
    m_sprimekey[1] = sprimemax;
}
