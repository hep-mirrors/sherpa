#include "BEAM/Main/Collider_Weight.H"

#include "ATOOLS/Org/Exception.H"

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

  m_rejection = ATOOLS::Settings::GetMainSettings()["BEAM_OVERLAP_REJECTION"]
                    .SetDefault(0)
                    .Get<int>();
  if (m_rejection == 0) return;

  // The rejection needs an impact parameter for each beam. Impact-parameter
  // integration variables are registered (in Beam_Channels) only for beams
  // carrying an EPA spectrum, so require at least one such beam
  if (p_beams[0]->Type() != beamspectrum::EPA &&
      p_beams[1]->Type() != beamspectrum::EPA)
    THROW(fatal_error,
          "BEAM_OVERLAP_REJECTION requires at least one EPA/Pomeron/Reggeon "
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

Collider_Weight::~Collider_Weight() = default;

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
  return m_weight;
}

double Collider_Weight::OverlapWeight()
{
  double b1(p_beams[0]->ImpactParameter()), b2(p_beams[1]->ImpactParameter());
  double b =
      std::sqrt(b1 * b1 + b2 * b2 - 2 * b1 * b2 * std::cos(2 * M_PI * m_eran));
  return (*p_rejector)(b);
}
