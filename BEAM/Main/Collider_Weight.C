#include "BEAM/Main/Collider_Weight.H"

#include "ATOOLS/Math/Random.H"

using namespace BEAM;

Collider_Weight::Collider_Weight(Kinematics_Base* kinematics)
    : Weight_Base(kinematics), m_mode(collidermode::unknown),
      p_rejector(nullptr)
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
                        .SetDefault(false)
                        .Get<int>();
  if (m_rejection == 0) return;
  if (m_rejection == 1)
    p_rejector = new Radius_Rejection(p_beams[0]->Beam(), p_beams[1]->Beam());
  else if (m_rejection > 1 && p_beams[0]->Beam().Kfcode() == kf_p_plus &&
           p_beams[1]->Beam().Kfcode() == kf_p_plus)
    p_rejector =
            new Proton_Proton_Rejection(p_beams[0]->Beam(), p_beams[1]->Beam());
  else if (m_rejection > 1 && (p_beams[0]->Beam().Kfcode() == kf_p_plus &&
                               p_beams[1]->Beam().IsIon()) ||
           (p_beams[0]->Beam().IsIon() &&
            p_beams[1]->Beam().Kfcode() == kf_p_plus))
    p_rejector = new Proton_Nucleon_Rejection(p_beams[0]->Beam(),
                                              p_beams[1]->Beam());
  else if (m_rejection > 1 && p_beams[0]->Beam().IsIon() &&
           p_beams[1]->Beam().IsIon())
    p_rejector = new Nucleon_Nucleon_Rejection(p_beams[0]->Beam(),
                                               p_beams[1]->Beam());
  else
    msg_Error() << METHOD << ATOOLS::om::red
                << ": Could not find appropriate beam radius rejection model, "
                   "will not do any rejection.\n"
                << ATOOLS::om::reset;
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
  if (m_rejection > 0 && Reject()) return 0.;
  m_weight = p_beams[0]->Weight() * p_beams[1]->Weight();
  return m_weight;
}

bool Collider_Weight::Reject()
{
  double b1(p_beams[0]->ImpactParameter()), b2(p_beams[1]->ImpactParameter());
  double b = std::sqrt(b1 * b1 + b2 * b2 -
                       2 * b1 * b2 * std::cos(2 * M_PI * ATOOLS::ran->Get()));
  return (*p_rejector)(b);
}
