#include "BEAM/Main/Transverse_Profile.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "BEAM/Main/Woods_Saxon.H"

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
  case Transverse_Profile_Type::WoodsSaxon: return os << "WoodsSaxon";
  }
  return os << "UnknownTransverseProfile";
}

std::istream& BEAM::operator>>(std::istream& is, Transverse_Profile_Type& type)
{
  std::string token;
  if (!(is >> token)) return is;
  if (token == "Gaussian") type = Transverse_Profile_Type::Gaussian;
  else if (token == "UniformSphere") type = Transverse_Profile_Type::UniformSphere;
  else if (token == "WoodsSaxon") type = Transverse_Profile_Type::WoodsSaxon;
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
    THROW(fatal_error, "A Woods-Saxon profile is constructed from its density.");
  }
  if (!std::isfinite(m_scale))
    THROW(fatal_error, "Transverse profile radius is too large.");
}

Transverse_Profile::Transverse_Profile(const Woods_Saxon& density)
    : m_type(Transverse_Profile_Type::WoodsSaxon),
      m_rms_radius(density.TransverseRMSRadius()), m_scale(1.)
{
  // In t = s^2 the radial CDF has the slope dF/dt = pi T(s), so F is the
  // running integral of the thickness T. Tabulated on a uniform grid in t, and
  // interpolated linearly in t, F is accurate to about 1e-6. Beyond
  // s = R + 20 d, 1 - F is below 1e-9.
  constexpr size_t nodes = 4000;
  const double smax = density.R() + 20. * density.d();
  ATOOLS::OneDim_Table thickness(ATOOLS::axis(nodes, 0., ATOOLS::sqr(smax)));
  for (size_t i = 0; i <= nodes; ++i)
    thickness.Fill(i, density.Thickness(std::sqrt(thickness.GetAxis().x(i))));
  double norm; // 1/pi, up to the truncation and the integration error
  auto cdf = thickness.Cumulative(0., norm);
  cdf->Rescale(1. / norm);
  p_cdf = std::move(cdf);
}

double Transverse_Profile::SampleRadius(double u) const
{
  if (!std::isfinite(u) || u < 0. || u >= 1.)
    THROW(fatal_error, "Transverse profile coordinate must lie in [0,1).");
  if (p_cdf) return std::sqrt(p_cdf->Inverse(u));
  if (m_scale == 0. || u == 0.) return 0.;
  const double log_survival = std::log1p(-u);
  if (m_type == Transverse_Profile_Type::Gaussian)
    return m_scale * std::sqrt(-log_survival);
  // expm1 avoids cancellation in 1 - (1-u)^(2/3) near u = 0.
  return m_scale * std::sqrt(-std::expm1((2. / 3.) * log_survival));
}
