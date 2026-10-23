#include "BEAM/Main/Collider_Weight.H"

#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "BEAM/Spectra/Gaussian.H"

#include <cmath>
#include <functional>

using namespace BEAM;

Collider_Weight::Collider_Weight(Kinematics_Base* kinematics)
    : Weight_Base(kinematics), m_mode(collidermode::unknown),
      p_rejector(nullptr), m_eran(0.)
{
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

  /*
    BEAM_SPREAD_CORRELATION rho: the two Gaussian beam energies follow the
    correlated two-dimensional Gaussian
  */
  m_rho = ATOOLS::Settings::GetMainSettings()["BEAM_SPREAD_CORRELATION"]
              .SetDefault(0.0)
              .Get<double>();
  if (m_rho != 0.) {
    for (int i = 0; i < 2; ++i) p_gauss[i] = dynamic_cast<const Gaussian*>(p_beams[i]);
    if (!p_gauss[0] || !p_gauss[1])
      THROW(fatal_error, "BEAM_SPREAD_CORRELATION needs BEAM_SPECTRA: "
                         "[Gaussian, Gaussian].");
    if (!(std::fabs(m_rho) < 1.))
      THROW(fatal_error, "BEAM_SPREAD_CORRELATION must lie in (-1, 1).");
    const double n1(p_gauss[0]->NSigma()), n2(p_gauss[1]->NSigma());
    m_rhonorm = std::erf(n1 / M_SQRT2) * std::erf(n2 / M_SQRT2)
                / BoxProbability(m_rho, n1, n2);
    msg_Info() << "Correlated Gaussian beam energy spread: rho = " << m_rho
               << ", box renormalisation " << m_rhonorm << ".\n";
  }

  m_rejection = ATOOLS::Settings::GetMainSettings()["BEAM_OVERLAP_REJECTION"]
                    .SetDefault(0)
                    .Get<int>();
  if (m_rejection <= 0) return;

  // The rejection needs an impact parameter for each beam. Impact-parameter
  // integration variables are registered (in Beam_Channels) only for beams
  // carrying an EPA spectrum, so require at least one such beam
  if (p_beams[0]->Type() != beamspectrum::EPA &&
      p_beams[1]->Type() != beamspectrum::EPA)
    THROW(fatal_error,
          "BEAM_OVERLAP_REJECTION requires at least one EPA "
          "beam to define an impact parameter.");

  const ATOOLS::Flavour& b0 = p_beams[0]->Beam();
  const ATOOLS::Flavour& b1 = p_beams[1]->Beam();
  if (m_rejection == 1)
    p_rejector = new Radius_Rejection(b0, b1);
  else if (b0.Kfcode() == kf_p_plus && b1.Kfcode() == kf_p_plus)
    p_rejector = new Proton_Proton_Rejection(
        b0, b1, (p_beams[0]->InMomentum() + p_beams[1]->InMomentum()).Abs2());
  else if ((b0.Kfcode() == kf_p_plus && b1.IsIon()) ||
           (b0.IsIon() && b1.Kfcode() == kf_p_plus)) {
    // fail fast at setup rather than aborting at the first weight evaluation
    THROW(not_implemented,
          "Proton-nucleon beam overlap rejection is not implemented.");
  } else if (b0.IsIon() && b1.IsIon()) {
    THROW(not_implemented,
          "Nucleon-nucleon beam overlap rejection is not implemented.");
  }

  if (p_rejector == nullptr)
    THROW(fatal_error,
          "BEAM_OVERLAP_REJECTION requested but no rejection model matches the "
          "beam combination.");
}

Collider_Weight::~Collider_Weight() { delete p_rejector; }

void Collider_Weight::AssignKeys(ATOOLS::Integration_Info* const info)
{
  m_sprimekey.Assign(m_keyid + std::string("s'"), 5, 0, info);
  m_ykey.Assign(m_keyid + std::string("y"), 3, 0, info);
  // Convention for m_xkey:
  // [x_{min,beam0}, x_{min,beam1}, x_{max,beam0}, x_{max,beam1}, x_{val,beam0},
  // x_{val,beam1}]. The limits, i.e. index 0,1,2,3 are saved as log(x), the
  // values are saved linearly.
  m_xkey.Assign(m_keyid + std::string("x"), 6, 0, info);
}

bool Collider_Weight::Calculate(const double& scale)
{
  m_weight = 0.;
  return (p_beams[0]->CalculateWeight(m_xkey[4], scale) &&
          p_beams[1]->CalculateWeight(m_xkey[5], scale));
}

double Collider_Weight::operator()()
{
  double overlap_weight(1.);
  if (m_rejection > 0) overlap_weight *= OverlapWeight();
  m_weight = p_beams[0]->Weight() * p_beams[1]->Weight() * overlap_weight;
  if (m_rho != 0.) m_weight *= CorrelationFactor();
  return m_weight;
}

double Collider_Weight::CorrelationFactor()
{
  const double d1((m_xkey[4] - p_gauss[0]->X0()) / p_gauss[0]->SigmaX());
  const double d2((m_xkey[5] - p_gauss[1]->X0()) / p_gauss[1]->SigmaX());
  return CorrelationWeight(d1, d2, m_rho, m_rhonorm);
}

double Collider_Weight::CorrelationWeight(double d1, double d2, double rho,
                                          double rhonorm)
{
  const double omr2(1. - rho * rho);
  return rhonorm / std::sqrt(omr2) *
         std::exp(-(rho * rho * (d1 * d1 + d2 * d2) - 2. * rho * d1 * d2)
                  / (2. * omr2));
}

/*
  Probability of the standard correlated Gaussian in |d1| < n1, |d2| < n2:
  the d2 integral is done analytically (normal CDF), the d1 integral with
  ATOOLS::Gauss_Integrator (Gauss-Legendre) to a relative precision of 1e-12.
*/
double Collider_Weight::BoxProbability(double rho, double n1, double n2)
{
  const double s(std::sqrt(1. - rho * rho));
  auto Phi = [](double z) { return 0.5 * std::erfc(-z / M_SQRT2); };
  const std::function<double(double)> f = [&](double d1) {
    return std::exp(-0.5 * d1 * d1) / std::sqrt(2. * M_PI) *
           (Phi((n2 - rho * d1) / s) - Phi((-n2 - rho * d1) / s));
  };
  ATOOLS::Lambda_Functor functor(&f);
  ATOOLS::Gauss_Integrator integrator(&functor);
  return integrator.Integrate(-n1, n1, 1.e-12);
}

double Collider_Weight::OverlapWeight()
{
  double b1(p_beams[0]->ImpactParameter()), b2(p_beams[1]->ImpactParameter());
  double b =
      std::sqrt(b1 * b1 + b2 * b2 - 2 * b1 * b2 * std::cos(2 * M_PI * m_eran));
  return (*p_rejector)(b);
}
