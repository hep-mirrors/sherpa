#include "BEAM/Main/Transverse_Profile.H"

#include "ATOOLS/Org/Exception.H"

#include <cmath>
#include <istream>
#include <ostream>
#include <string>

using namespace BEAM;

std::ostream& BEAM::operator<<(std::ostream& os, Transverse_Profile_Type type)
{
  switch (type) {
  case Transverse_Profile_Type::Gaussian: return os << "Gaussian";
  case Transverse_Profile_Type::UniformSphere: return os << "UniformSphere";
  }
  return os << "UnknownTransverseProfile";
}

std::istream& BEAM::operator>>(std::istream& is, Transverse_Profile_Type& type)
{
  std::string token;
  if (!(is >> token)) return is;
  if (token == "Gaussian") type = Transverse_Profile_Type::Gaussian;
  else if (token == "UniformSphere") type = Transverse_Profile_Type::UniformSphere;
  else is.setstate(std::ios::failbit);
  return is;
}

Transverse_Profile::Transverse_Profile(Transverse_Profile_Type type,
                                     double rms_radius)
    : m_type(type), m_rms_radius(rms_radius), m_scale(rms_radius)
{
  if (!std::isfinite(rms_radius) || rms_radius < 0.)
    THROW(fatal_error, "Transverse profile RMS radius must be finite and non-negative.");
  switch (type) {
  case Transverse_Profile_Type::Gaussian: break;
  case Transverse_Profile_Type::UniformSphere:
    m_scale *= std::sqrt(2.5); // outer radius; <s^2> = 2 R^2 / 5
    break;
  default:
    THROW(fatal_error, "Unknown transverse profile type.");
  }
  if (!std::isfinite(m_scale))
    THROW(fatal_error, "Transverse profile radius is too large.");
}

double Transverse_Profile::SampleRadius(double u) const
{
  if (!std::isfinite(u) || u < 0. || u >= 1.)
    THROW(fatal_error, "Transverse profile coordinate must lie in [0,1).");
  if (m_scale == 0. || u == 0.) return 0.;

  const double log_survival = std::log1p(-u);
  if (m_type == Transverse_Profile_Type::Gaussian)
    return m_scale * std::sqrt(-log_survival);
  // expm1 avoids cancellation in 1 - (1-u)^(2/3) near u = 0.
  return m_scale * std::sqrt(-std::expm1((2. / 3.) * log_survival));
}
