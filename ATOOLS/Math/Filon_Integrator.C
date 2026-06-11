#include "ATOOLS/Math/Filon_Integrator.H"

#include "ATOOLS/Org/Exception.H"

#include <cmath>
#include <utility>

using namespace ATOOLS;

Filon_Integrator::Filon_Integrator(Function_Base* f, double a, double b,
                                   size_t npanels)
    : m_a(a), m_b(b), m_npanels(npanels)
{
  if (npanels < 1) THROW(fatal_error, "Filon_Integrator needs >= 1 panel.");
  if (b <= a) THROW(fatal_error, "Filon_Integrator needs b > a.");
  const size_t nnodes = 2 * npanels + 1;
  m_h = (b - a) / double(nnodes - 1);
  m_f.resize(nnodes);
  for (size_t j = 0; j < nnodes; ++j) m_f[j] = (*f)(a + double(j) * m_h);
}

Filon_Integrator::Filon_Integrator(std::vector<double> fvals, double a, double b)
    : m_a(a), m_b(b), m_f(std::move(fvals))
{
  if (m_f.size() < 3 || m_f.size() % 2 == 0)
    THROW(fatal_error, "Filon_Integrator needs an odd number (>= 3) of nodes.");
  if (b <= a) THROW(fatal_error, "Filon_Integrator needs b > a.");
  m_npanels = (m_f.size() - 1) / 2;
  m_h = (b - a) / double(m_f.size() - 1);
}

void Filon_Integrator::Coefficients(double theta, double& alpha, double& beta,
                                    double& gamma)
{
  if (std::abs(theta) < 0.1) {
    // Small-theta Taylor expansions (A&S 25.4.47) avoid the catastrophic
    // cancellation of the closed form (alpha is an O(theta^3) difference of
    // O(theta) terms).
    const double t2 = theta * theta;
    alpha = theta * t2 * (2. / 45. - t2 * (2. / 315. - t2 * (2. / 4725.)));
    beta = 2. / 3. + t2 * (2. / 15. - t2 * (4. / 105. - t2 * (2. / 567.)));
    gamma = 4. / 3. - t2 * (2. / 15. - t2 * (1. / 210. - t2 * (1. / 11340.)));
  } else {
    const double s = std::sin(theta), c = std::cos(theta);
    const double t3 = theta * theta * theta;
    alpha = (theta * theta + theta * s * c - 2. * s * s) / t3;
    beta = 2. * (theta * (1. + c * c) - 2. * s * c) / t3;
    gamma = 4. * (s - theta * c) / t3;
  }
}

double Filon_Integrator::SineTransform(double k) const
{
  if (k == 0.) return 0.;
  double alpha, beta, gamma;
  Coefficients(k * m_h, alpha, beta, gamma);

  const size_t last = m_f.size() - 1;
  // Even-index nodes, endpoints counted with weight 1/2.
  double s_even = -0.5 * (m_f[0] * std::sin(k * m_a) +
                          m_f[last] * std::sin(k * m_b));
  for (size_t j = 0; j <= last; j += 2)
    s_even += m_f[j] * std::sin(k * (m_a + double(j) * m_h));
  // Odd-index nodes.
  double s_odd = 0.;
  for (size_t j = 1; j < last; j += 2)
    s_odd += m_f[j] * std::sin(k * (m_a + double(j) * m_h));

  return m_h * (alpha * (m_f[0] * std::cos(k * m_a) -
                         m_f[last] * std::cos(k * m_b)) +
                beta * s_even + gamma * s_odd);
}

double Filon_Integrator::CosineTransform(double k) const
{
  const size_t last = m_f.size() - 1;
  if (k == 0.) {
    // Plain \int_a^b f dr via composite Simpson on the same grid.
    double s = m_f[0] + m_f[last];
    for (size_t j = 1; j < last; j += 2) s += 4. * m_f[j];
    for (size_t j = 2; j < last; j += 2) s += 2. * m_f[j];
    return m_h / 3. * s;
  }
  double alpha, beta, gamma;
  Coefficients(k * m_h, alpha, beta, gamma);

  double c_even = -0.5 * (m_f[0] * std::cos(k * m_a) +
                          m_f[last] * std::cos(k * m_b));
  for (size_t j = 0; j <= last; j += 2)
    c_even += m_f[j] * std::cos(k * (m_a + double(j) * m_h));
  double c_odd = 0.;
  for (size_t j = 1; j < last; j += 2)
    c_odd += m_f[j] * std::cos(k * (m_a + double(j) * m_h));

  return m_h * (alpha * (m_f[last] * std::sin(k * m_b) -
                         m_f[0] * std::sin(k * m_a)) +
                beta * c_even + gamma * c_odd);
}
