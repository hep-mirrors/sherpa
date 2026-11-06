#include "BEAM/Main/Transverse_Profile.H"

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

#include <algorithm>
#include <cmath>
#include <functional>
#include <istream>
#include <ostream>
#include <string>

using namespace BEAM;

namespace {

  constexpr double skin_range = 40.;
  constexpr double quadrature_tolerance = 1.e-10;
  constexpr double cdf_tolerance = 1.e-8;

  // Integral dz z^power g(hypot(s,z)). Below R0-40a, g=1 to floating-point
  // accuracy; beyond R0+40a the omitted probability is below 4e-15 even for
  // R0=0. Splitting off the core also resolves arbitrarily thin surfaces.
  double WoodsSaxonIntegral(unsigned int power, double radius, double skin,
                            double s = 0.)
  {
    const double outer = radius + skin_range * skin;
    if (s >= outer) return 0.;
    const double inner = std::max(0., radius - skin_range * skin);
    const double zmin = s < inner ? std::sqrt((inner - s) * (inner + s)) : 0.;
    const double zmax = std::sqrt((outer - s) * (outer + s));
    const std::function<double(double)> integrand = [&](double z) {
      return std::pow(z, power) /
             (1. + std::exp((std::hypot(s, z) - radius) / skin));
    };
    ATOOLS::Lambda_Functor function(&integrand);
    ATOOLS::Gauss_Integrator integral(&function);
    return std::pow(zmin, power + 1) / (power + 1) +
           integral.Integrate(zmin, zmax, quadrature_tolerance);
  }

  // Transverse <s^2>. Rescale before taking moments so the R0=0 bracket
  // remains well-conditioned even when skin/RMS is extremely small.
  double WoodsSaxonMeanSquare(double radius, double skin)
  {
    const double scale = std::max(radius, skin);
    return (2. / 3.) * scale * scale *
        WoodsSaxonIntegral(4, radius / scale, skin / scale) /
        WoodsSaxonIntegral(2, radius / scale, skin / scale);
  }

} // namespace

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
                                     double rms_radius, double skin_depth)
    : m_type(type), m_rms_radius(rms_radius), m_scale(rms_radius),
      m_skin_depth(type == Transverse_Profile_Type::WoodsSaxon ? skin_depth : 0.)
{
  if (!std::isfinite(rms_radius) || rms_radius < 0.)
    THROW(fatal_error, "Transverse profile RMS radius must be finite and non-negative.");
  switch (type) {
  case Transverse_Profile_Type::Gaussian: break;
  case Transverse_Profile_Type::UniformSphere:
    m_scale *= std::sqrt(2.5); // outer radius; <s^2> = 2 R^2 / 5
    break;
  case Transverse_Profile_Type::WoodsSaxon:
    if (!std::isfinite(skin_depth) || skin_depth < 0.)
      THROW(fatal_error, "Woods-Saxon skin depth must be finite and non-negative.");
    if (skin_depth == 0.) m_scale *= std::sqrt(2.5);
    else if (rms_radius > 0.) BuildWoodsSaxon();
    break;
  default:
    THROW(fatal_error, "Unknown transverse profile type.");
  }
  if (!std::isfinite(m_scale))
    THROW(fatal_error, "Transverse profile radius is too large.");
}

void Transverse_Profile::BuildWoodsSaxon()
{
  // Work in units of the requested transverse RMS radius. Determine R0 so
  // that the projected density, not an asymptotic formula, has <s^2>=1.
  const double skin = m_skin_depth / m_rms_radius;
  if (skin == 0.) {
    m_scale *= std::sqrt(2.5); // diffuseness below floating-point resolution
    return;
  }
  if (!std::isfinite(skin) || skin >= 1. ||
      WoodsSaxonMeanSquare(0., skin) > 1.)
    THROW(fatal_error, "Woods-Saxon RMS radius is too small for this skin depth.");
  double lo = 0., hi = std::sqrt(2.5);
  for (unsigned int i = 0; i < 40; ++i) {
    const double mid = .5 * (lo + hi);
    if (WoodsSaxonMeanSquare(mid, skin) < 1.) lo = mid;
    else hi = mid;
  }
  const double radius = .5 * (lo + hi);
  const double norm = WoodsSaxonIntegral(2, radius, skin);
  if (!(norm > 0.) || !std::isfinite(norm))
    THROW(fatal_error, "Could not normalise the Woods-Saxon profile.");

  // 1-F(s) = integral dz z^2 g(hypot(s,z)) / integral dr r^2 g(r).
  // This gives the projected radial CDF without a nested thickness integral.
  const auto cdf = [&](double s2) {
    return std::clamp(1. - WoodsSaxonIntegral(2, radius, skin, std::sqrt(s2)) /
                              norm, 0., 1.);
  };
  m_inverse_cdf.emplace_back(0., 0.);
  const auto refine = [&](auto&& self, double t0, double u0,
                          double t1, double u1, unsigned int depth) -> void {
    const double tm = .5 * (t0 + t1);
    const double um = std::clamp(cdf(tm), u0, u1);
    // F is concave in s^2: its maximum linear-interpolation error is at most
    // twice the midpoint error. Interpolating s^2 also preserves F~s^2 at 0.
    if (std::abs(um - .5 * (u0 + u1)) > .5 * cdf_tolerance) {
      if (depth == 24 || m_inverse_cdf.size() >= 65536)
        THROW(fatal_error, "Woods-Saxon CDF table did not converge.");
      self(self, t0, u0, tm, um, depth + 1);
      self(self, tm, um, t1, u1, depth + 1);
    } else if (u1 > m_inverse_cdf.back().first) {
      m_inverse_cdf.emplace_back(u1, t1);
    }
  };
  const double outer = radius + skin_range * skin;
  refine(refine, 0., 0., outer * outer, 1., 0);
}

double Transverse_Profile::SampleRadius(double u) const
{
  if (!std::isfinite(u) || u < 0. || u >= 1.)
    THROW(fatal_error, "Transverse profile coordinate must lie in [0,1).");
  if (m_scale == 0. || u == 0.) return 0.;
  if (!m_inverse_cdf.empty()) {
    const auto upper = std::upper_bound(
        m_inverse_cdf.begin(), m_inverse_cdf.end(), u,
        [](double value, const auto& node) { return value < node.first; });
    const auto lower = upper - 1;
    const double fraction = (u - lower->first) / (upper->first - lower->first);
    return m_scale * std::sqrt(lower->second +
                               fraction * (upper->second - lower->second));
  }

  const double log_survival = std::log1p(-u);
  if (m_type == Transverse_Profile_Type::Gaussian)
    return m_scale * std::sqrt(-log_survival);
  // expm1 avoids cancellation in 1 - (1-u)^(2/3) near u = 0.
  return m_scale * std::sqrt(-std::expm1((2. / 3.) * log_survival));
}
