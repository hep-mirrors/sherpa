#include "ATOOLS/Math/Bessel_Integrator.H"
#include <iomanip>

using namespace ATOOLS;

Bessel_Integrator::Bessel_Integrator(ATOOLS::Function_Base* f,
                                     const size_t& order)
    : m_kernel(f, order), m_order(order), m_maxbins(300), m_depth(10),
      m_iterator(1)
{
  FixBins(false);
  m_F.resize(m_maxbins, 0.);
  m_Psi.resize(m_maxbins, 0.);
  m_M.resize(m_depth);
  m_N.resize(m_depth);
  for (size_t i = 0; i < m_depth; i++) {
    m_M[i].resize(m_maxbins, 0.);
    m_N[i].resize(m_maxbins, 0.);
  }
  m_gauss = Gauss_Integrator(&m_kernel);
}

double Bessel_Integrator::operator()(double tolerance, bool output)
{
  size_t stability_counter = 0;

  if (output)
    msg_Out() << "=== " << METHOD << " start filling the supports, "
              << "depth = " << m_depth << ", maxbins = " << m_maxbins << "\n";

  double F = 0.;
  double xmin, xmax;
  for (size_t i = 1; i < m_maxbins + 1; i++) {
    if (m_iterator == 0) {
      xmin = m_zeroes[i - 1];
      xmax = m_zeroes[i];
    } else {
      xmin = m_extrema[i - 1];
      xmax = m_extrema[i];
    }
    F += m_Psi[i - 1] = m_gauss.Integrate(xmin, xmax, tolerance);
    m_F[i - 1] = F;
    m_M[0][i - 1] = m_F[i - 1] / m_Psi[i - 1];
    m_N[0][i - 1] = 1. / m_Psi[i - 1];
  if (output) {
      msg_Out() << "  bin i = " << std::setw(2) << i << ": "
                << "Psi[" << std::setw(8) << xmin << ", " << std::setw(8)
                << xmax << "] = " << std::setw(12) << std::setprecision(6)
                << m_Psi[i - 1] << ", "
                << "F[" << std::setw(8) << xmin << ", " << std::setw(8) << xmax
                << "] = " << std::setw(12) << std::setprecision(6) << m_F[i - 1]
                << " --> "
                << "M_{-1}^(" << std::setw(2) << i << ") = " << std::setw(12)
                << std::setprecision(6) << m_M[0][i - 1] << ", and "
                << "N_{-1}^(" << std::setw(2) << i << ") = " << std::setw(12)
                << std::setprecision(6) << m_N[0][i - 1] << "\n";
    }
  }

  // adapted the naming from eq. 1.12 and 1.13 in DOI:10.2307/2008589
  for (size_t p = 1; p < m_depth; p++) {
    for (size_t s = 1; s < m_maxbins + 1; s++) {
      if (m_iterator == 0) {
        xmin = m_zeroes[s];
        xmax = s + 1 + p > m_maxbins ? m_zeroes[s] + (p + 1) * M_PI
                                     : m_zeroes[s + 1 + p];
      } else {
        xmin = m_extrema[s];
        xmax = s + 1 + p > m_maxbins ? m_extrema[s] + (p + 1) * M_PI
                                     : m_extrema[s + 1 + p];
      }
      m_M[p][s - 1] =
          (m_M[p - 1][s - 1] - m_M[p - 1][s]) * xmin * xmax / (xmax - xmin);
      m_N[p][s - 1] =
          (m_N[p - 1][s - 1] - m_N[p - 1][s]) * xmin * xmax / (xmax - xmin);
    }

    // Convergence check
    if (p >= 3 && std::abs(m_N[p][0]) > 1.e-30 && std::abs(m_N[p-1][0]) > 1.e-30) {
      double current_result = m_M[p][0] / m_N[p][0];
      double prev_result = m_M[p - 1][0] / m_N[p - 1][0];
      double rel_error = std::abs(current_result - prev_result)
                         / (std::abs(current_result) + 1.0e-20);

      if (rel_error < tolerance) {
        stability_counter++;
      } else {
        stability_counter = 0;
      }

      if (stability_counter >= 2) {
        if (output) {
          msg_Out() << "=== " << METHOD << ": Converged at depth " << p
                    << ", Result = " << std::setprecision(8) << m_M[p - 1][0] / m_N[p - 1][0] << "\n";
        }
        return m_M[p - 1][0] / m_N[p - 1][0];
      }
    }
  }

  // Fallback
  if (output) {
    msg_Out() << "=== " << METHOD << ": Did NOT converge within depth=" << m_depth
              << "\n";
    msg_Out() << "    Will fall back to cumulatd sum: " << std::setprecision(16)
              << m_F.back() << "\n";
  }
  return m_F.back();
}

void Bessel_Integrator::FixBins(const bool& output)
{
  m_zeroes.clear();
  m_zeroes.reserve(m_maxbins + 1);
  m_extrema.resize(m_maxbins + 1);
  // We fix the first value at zero anyway, as the integration
  // always starts at 0 for our purposes
  switch (m_order) {
  case 0:
    m_zeroes.assign({0.,
                     2.404825557695773,
                     5.520078110286311,
                     8.653727912911012,
                     11.79153443901428,
                     14.93091770848779,
                     18.07106396791092,
                     21.21163662987926,
                     24.35247153074930,
                     27.49347913204025,
                     30.63460646843198,
                     33.77582021357357,
                     36.91709835366404,
                     40.05842576462824,
                     43.19979171317673,
                     46.34118837166181,
                     49.48260989739782,
                     52.62405184111500,
                     55.76551075501998,
                     58.90698392608094,
                     62.04846919022717,
                     65.18996480020686,
                     68.33146932985680,
                     71.47298160359373,
                     74.61450064370184,
                     77.75602563038806,
                     80.89755587113763,
                     84.03909077693819,
                     87.18062984364115,
                     90.32217263721048,
                     93.46371878194477});
    break;
  case 1:
    m_zeroes.assign({0.,
                     3.831705970207515,
                     7.015586669815619,
                     10.17346813506272,
                     13.32369193631422,
                     16.47063005087763,
                     19.61585851046824,
                     22.76008438059277,
                     25.90367208761838,
                     29.04682853491686,
                     32.18967991097440,
                     35.33230755008387,
                     38.47476623477162,
                     41.61709421281445,
                     44.75931899765282,
                     47.90146088718545,
                     51.04353518357151,
                     54.18555364106132,
                     57.32752543790101,
                     60.46945784534749,
                     63.61135669848123,
                     66.75322673409849,
                     69.89507183749577,
                     73.03689522557383,
                     76.17869958464146,
                     79.32048717547630,
                     82.46225991437356,
                     85.60401943635023,
                     88.74576714492631,
                     91.88750425169499,
                     95.02923180804470});
    break;
    /*
     * these are not needed as in fact the higher order Bessel functions are
     * not implemented in our code. Are left here in case one day we need to
     * extend. case 2:
     *   m_zeroes.assign({0., 5.1356, 8.4172, 11.6198, 14.7960, 17.9598});
     *   break;
     * default:
     *     m_zeroes.assign({0.,
     *                      (m_order + 1.8557571 * pow(m_order, 1. / 3.) +
     *                    1.033150 * pow(m_order, -1. / 3.) - 0.00397 /
     * m_order - 0.0908 * pow(m_order, -5. / 3.) + 0.043 * pow(m_order, -7.
     * / 3.)), (m_order + 3.2446076 * pow(m_order, 1. / 3.) + 3.158244 *
     * pow(m_order, -1. / 3.) - 0.08331 / m_order - 0.8437 * pow(m_order, -5.
     * / 3.) + 0.864 * pow(m_order, -7. / 3.))}); break;
     */
  }
  // Fill in the remaining zeroes
  for (int i = m_zeroes.size(); i < m_maxbins + 1; ++i)
    m_zeroes.push_back(m_zeroes[i - 1] + M_PI);
  // approximate the extrema as midpoints between zeroes; again, the
  // integration starts at 0, so we fix the first value there
  m_extrema[0] = 0.;
  for (int i = 1; i < m_maxbins + 1; ++i) {
    m_extrema[i] = (m_zeroes[i] + m_zeroes[i - 1]) / 2.;
  }
  if (output) {
    msg_Out() << "=== " << METHOD << " yields the first " << m_maxbins << " "
              << "zeroes and extrema of Bessel J_{" << m_order << "}:\n";
    for (size_t i = 0; i < m_maxbins + 1; i++) {
      msg_Out() << "  x_{" << m_order << ", " << std::setw(2) << i
                << "} = " << std::setw(12) << std::setprecision(8)
                << m_zeroes[i] << " --> "
                << "y_{" << m_order << ", " << std::setw(2) << i
                << "} = " << std::setw(12) << std::setprecision(8)
                << m_extrema[i] << ".\n";
    }
  }
}
