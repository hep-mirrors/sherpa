#include "ATOOLS/Math/Bessel_Integrator.H"
#include <iomanip>

using namespace ATOOLS;

Bessel_Integrator::Bessel_Integrator(ATOOLS::Function_Base* f,
                                     const size_t&          order)
    : m_kernel(f, order), m_order(order), m_maxbins(20), m_depth(10),
      m_maxdepth(10), m_iterator(1)
{
  FixBins(false);
  m_F.resize(m_maxbins + 1, 0.);
  m_Psi.resize(m_maxbins + 1, 0.);
  m_M.resize(m_depth);
  m_N.resize(m_depth);
  for (size_t i = 0; i < m_depth; i++) {
    m_M[i].resize(m_maxbins + 1, 0.);
    m_N[i].resize(m_maxbins + 1, 0.);
  }
}

double Bessel_Integrator::operator()()
{
  if (FillBins(false)) return (m_M[m_maxdepth - 1][0] / m_N[m_maxdepth - 1][0]);
  return 0.;
}

bool Bessel_Integrator::FillBins(const bool& output)
{
  double tolerance = 1.e-16;

  if (output)
    msg_Out() << "=== " << METHOD << " start filling the supports, "
              << "max depth = " << m_maxdepth << ":\n";
  Gauss_Integrator gauss(&m_kernel);
  double           F = 0.;
  double           xmin, xmax;
  for (size_t i = 1; i < m_maxbins + 1; i++) {
    if (m_iterator == 0) {
      xmin = m_zeroes[i - 1];
      xmax = m_zeroes[i];
    } else if (m_iterator == 1) {
      xmin = m_extrema[i - 1];
      xmax = m_extrema[i];
    } else
      exit(1);
    F += m_Psi[i - 1] = gauss.Integrate(xmin, xmax, 1.e-6);
    if (dabs(m_Psi[i - 1]) < tolerance && i < m_maxdepth + 1) {
      if (output)
        msg_Out() << "   found a zero in integral over "
                  << "x in [" << xmin << ", " << xmax << "]: " << m_Psi[i - 1]
                  << " for "
                  << "i = " << i << " --> new max depth = " << (i - 1) << "\n";
      m_maxdepth = i == 1 ? 1 : i - 1;
      return false;
    }
    m_F[i - 1] = F;
    if (dabs(m_Psi[i - 1]) > tolerance) {
      m_M[0][i - 1] = m_F[i - 1] / m_Psi[i - 1];
      m_N[0][i - 1] = 1. / m_Psi[i - 1];
    } else {
      m_M[0][i - 1] = 0.;
      m_N[0][i - 1] = 0.;
    }
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
  for (size_t p = 1; p < m_maxdepth; p++) {
    for (size_t s = 1; s < m_maxbins + 1; s++) {
      if (m_iterator == 0) {
        xmin = m_zeroes[s];
        xmax = s + 1 + p > m_maxbins ? m_zeroes[s] + (p + 1) * M_PI
                                     : m_zeroes[s + 1 + p];
      } else if (m_iterator == 1) {
        xmin = m_extrema[s];
        xmax = s + 1 + p > m_maxbins ? m_extrema[s] + (p + 1) * M_PI
                                     : m_extrema[s + 1 + p];
      } else
        exit(1);
      m_M[p][s - 1] = (m_M[p - 1][s - 1] - m_M[p - 1][s]) * xmin * xmax / (xmax - xmin);
      m_N[p][s - 1] = (m_N[p - 1][s - 1] - m_N[p - 1][s]) * xmin * xmax / (xmax - xmin);
    }
    if (output && p != m_maxdepth - 1) {
      msg_Out() << "=== " << METHOD << "(depth = " << std::setw(2) << p
                << "): " << std::setw(12) << std::setprecision(6)
                << (m_M[p][0] / m_N[p][0]) << "\n";
    }
  }
  if (output)
    msg_Out() << "=== " << METHOD << ": Integral value("
              << "depth = " << std::setw(2) << m_maxdepth
              << "): " << std::setw(12) << std::setprecision(6)
              << (m_M[m_maxdepth - 1][0] / m_N[m_maxdepth - 1][0]) << " "
              << "for f(x) * J_{" << m_order << "}(x)\n";
  return true;
}

void Bessel_Integrator::FixBins(const bool& output)
{
  m_zeroes.resize(m_maxbins + 1);
  m_extrema.resize(m_maxbins + 1);
  // We fix the first value at zero anyway, as the integration
  // always starts at 0 for our purposes
  switch (m_order) {
    case 0:
      m_zeroes.assign({0., 2.404825, 5.520078, 8.653727, 11.791534, 14.930917,
                       18.071063, 21.211636, 24.352471, 27.493479, 30.634606});
      break;
    case 1:
      m_zeroes.assign({0., 3.831705, 7.015586, 10.173468, 13.323691, 16.470630,
                       19.615858, 22.760084, 25.903672, 29.046828, 32.189679});
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
  for (int i = 11; i < m_maxbins + 1; ++i)
    m_zeroes[i] = m_zeroes[i - 1] + M_PI;
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

double Bessel_Integrator::TestFunction::operator()(double x)
{
  switch (m_mode) {
    case 4: return (1. - exp(-x)) / (x * log(1. + sqrt(2.)));
    case 3: return x / (1. + x * x);
    case 2: return log(1. + x * x) / 2.;
    case 1: return x * x;
    case 0:
    default: break;
  }
  return 1;
}

void Bessel_Integrator::Test()
{
  Function_Base* f;
  for (size_t i = 0; i < 5; i++) {
    m_kernel.SetFunc(f = new Bessel_Integrator::TestFunction(i));
    m_kernel.SetOrder((i == 2 || i == 0) ? 1 : 0);
    FillBins(true);
    delete f;
  }
}
