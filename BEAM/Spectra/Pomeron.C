#include "BEAM/Spectra/Pomeron.H"
#include "ATOOLS/Math/Random.H"

using namespace BEAM;
using namespace ATOOLS;

Pomeron::Pomeron(const Flavour _beam, const double _energy, const double _pol,
         const int _dir)
    : Beam_Base(beamspectrum::Pomeron, _beam, _energy, _pol, _dir), m_A(0.)
{
  Settings &s = Settings::GetMainSettings();
  m_proton_mass = Flavour(kf_p_plus).Mass();
  m_Nbunches = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_pomeron);
  m_bunches[1] = m_beam;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0.,0.,0.,0.);
  m_on = true;

  std::vector<double> tMax{s["Pomeron"]["tMax"].GetVector<double>()};
  if (tMax.size() != 1 && tMax.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `Pomeron:tMax'.");
  m_tMax = (_dir > 0) ? tMax.front() : tMax.back();

  std::vector<double> xMax{s["Pomeron"]["xMax"].GetVector<double>()};
  if (xMax.size() != 1 && xMax.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `Pomeron:xMax'.");
  m_xMax = (_dir > 0) ? xMax.front() : xMax.back();

  std::vector<double> B{s["Pomeron"]["B"].GetVector<double>()};
  if (B.size() != 1 && B.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `Pomeron:B'.");
  m_B = (_dir > 0) ? B.front() : B.back();

  std::vector<double> alpha_i{s["Pomeron"]["Alpha_intercept"].GetVector<double>()};
  if (alpha_i.size() != 1 && alpha_i.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `Pomeron:Alpha_intercept'.");
  m_alpha_intercept = (_dir > 0) ? alpha_i.front() : alpha_i.back();

  std::vector<double> alpha_s{s["Pomeron"]["Alpha_slope"].GetVector<double>()};
  if (alpha_s.size() != 1 && alpha_s.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `Pomeron:Alpha_slope'.");
  m_alpha_slope = (_dir > 0) ? alpha_s.front() : alpha_s.back();

  FixNormalisation();
}

bool Pomeron::CalculateWeight(double x, double q2)
{
  m_x = x;
  m_Q2 = q2;
  if (x > 1. - m_proton_mass / 2. / m_energy || x > m_xMax) {
    m_weight = 0.;
    return true;
  }
  double tmax = Min(2. * m_energy * (m_energy - m_proton_mass), m_tMax);
  double tmin = sqr(m_proton_mass * x) / (1 - x);
  /*
   * In analogy to the EPA spectrum, we integrated the original weight from
   * Goharipour:2018yov over t \in [-t_max, -tmin]; note that tmax is positive
   */
  m_weight =
          (m_A * std::pow(x, 1. - 2. * m_alpha_intercept) *
           (std::exp(-m_B * tmin) * std::pow(x, 2. * m_alpha_slope * tmin) -
            std::exp(-m_B * tmax) * std::pow(x, 2. * m_alpha_slope * tmax))) /
          (m_B - 2. * m_alpha_slope * std::log(x));
  if (m_weight < 0.) m_weight = 0.;
  return true;
}

void Pomeron::SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i) {
  if (i==0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab-out;
  }
}

void Pomeron::FixPosition() {
  // Sample the position of the pomeron uniformly across the proton
  double radius = Flavour(kf_p_plus).Radius() * ran->Get();
  double phi = 2*M_PI*ran->Get();
  m_position = Vec4D(0., radius*std::cos(phi), radius*std::sin(phi), 0.);
}

void Pomeron::FixNormalisation()
{
  // according to Vadim Guzey, in Goharipour:2018yov they normalized the flux as
  // in hep-ex/0606004, after eq. 14
  double x    = 0.003;
  double tmax = Min(2. * m_energy * (m_energy - m_proton_mass), m_tMax);
  double tmin = sqr(m_proton_mass * x) / (1 - x);
  m_A         = (m_B - 2. * m_alpha_slope * std::log(x)) /
        (std::pow(x, 2. - 2. * m_alpha_intercept) *
         (std::exp(-m_B * tmin) * std::pow(x, 2. * m_alpha_slope * tmin) -
          std::exp(-m_B * tmax) * std::pow(x, 2. * m_alpha_slope * tmax)));
}

Beam_Base *Pomeron::Copy() { return new Pomeron(*this); }
