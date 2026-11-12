#include "BEAM/Main/Woods_Saxon.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/Flavour.H"

#include <cmath>

using namespace BEAM;
using ATOOLS::sqr;

// The closed forms for the normalisation and <r^2> are the leading terms of the
// Sommerfeld expansion. They neglect corrections of relative size exp(-R/d):
// about 3e-4 for carbon and 1e-8 for lead with the default R and d.

Woods_Saxon::Woods_Saxon(double R, double d) : m_R(R), m_d(d)
{
  if (!std::isfinite(R) || !std::isfinite(d) || !(R > 0.) || !(d > 0.))
    THROW(fatal_error, "Woods-Saxon radius and skin depth must be finite and > 0.");
  // \int d^3r / (1 + exp((r - R) / d)) = 4 pi / 3 (R^3 + pi^2 d^2 R)
  m_rho0 = 3. / (4. * M_PI * R * (R * R + sqr(M_PI * d)));
}

Woods_Saxon Woods_Saxon::FromSettings(const ATOOLS::Flavour& nucleus,
                                      size_t beam)
{
  // The default radius depends on the nucleus, which may differ between the
  // beams, so it cannot be a common default of the setting.
  auto s = ATOOLS::Settings::GetMainSettings()["EPA"];
  s.DeclareVectorSettingsWithEmptyDefault({"WoodsSaxon_R"});
  double R = 1.118 * std::cbrt(nucleus.GetMassNumber());
  if (!s["WoodsSaxon_R"].GetVector<double>().empty())
    R = s["WoodsSaxon_R"].GetTwoVector<double>()[beam];
  else if (!nucleus.IsIon())
    THROW(fatal_error, "Set EPA:WoodsSaxon_R for the Woods-Saxon density of " +
                           nucleus.IDName() + ".");
  const double d = s["WoodsSaxon_d"].SetDefault(0.54).GetTwoVector<double>()[beam];
  return Woods_Saxon(R / ATOOLS::rpa->hBar_c(), d / ATOOLS::rpa->hBar_c());
}

double Woods_Saxon::Density(double r) const
{
  return m_rho0 / (1. + std::exp((r - m_R) / m_d));
}

double Woods_Saxon::Thickness(double s) const
{
  // The trapezoidal rule on a uniform grid in z converges exponentially for
  // this smooth and decaying integrand: four nodes per skin depth give about
  // 1e-8. The density is cut at R + 40 d, where it has fallen by exp(-40).
  const double h = m_d / 4., rmax2 = sqr(m_R + 40. * m_d);
  if (s * s >= rmax2) return 0.;
  double sum = .5 * Density(s); // the integrand is even in z
  for (int k = 1; s * s + sqr(k * h) < rmax2; ++k)
    sum += Density(std::hypot(s, k * h));
  return 2. * h * sum;
}

double Woods_Saxon::TransverseRMSRadius() const
{
  // <r^2> = 3/5 (R^4 + 10/3 pi^2 d^2 R^2 + 7/3 pi^4 d^4) / (R^2 + pi^2 d^2)
  const double R2 = sqr(m_R), p2 = sqr(M_PI * m_d);
  return std::sqrt(.4 * (sqr(R2) + 10. / 3. * p2 * R2 + 7. / 3. * sqr(p2)) /
                   (R2 + p2));
}
