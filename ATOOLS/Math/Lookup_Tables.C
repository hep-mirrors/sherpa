#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Org/Exception.H"
#include <fstream>
#include <iostream>
#include <algorithm> // For std::lower_bound, std::min
#include <cmath>     // For std::fabs, std::sqrt, std::exp, etc.

using namespace ATOOLS;

axis::axis(size_t nbins, double xmin, double xmax,
           axis_mode::code mode)
    : m_nbins(nbins), m_xmin(xmin), m_xmax(xmax), m_mode(mode), m_xstep(0.0)
{
  if (m_nbins < 1)
    THROW(fatal_error, "Cannot create table with no bins.");
  if (xmin > xmax)
    THROW(fatal_error, "X_min must be smaller than X_max for table.");
  if (m_mode == axis_mode::linear) {
    m_xstep = (m_xmax - m_xmin) / double(m_nbins);
  } else if (m_mode == axis_mode::log) {
    if (m_xmin <= 0.0 || m_xmax <= 0.0)
      THROW(fatal_error, "Log axis requires positive min/max.");
    m_xstep = log(m_xmax / m_xmin) / double(m_nbins);
  }
}

double axis::x(size_t bin) const
{
  if (m_nbins == 0) return m_xmin;
  if (bin > m_nbins) {
    // This case should ideally not be reached by interpolation logic
    return m_xmax;
  }
  if (m_mode == axis_mode::linear) return m_xmin + double(bin) * m_xstep;
  // Assumes axis_mode::log
  return m_xmin * std::exp(m_xstep * double(bin));
}

size_t axis::bin(double x) const
{
  if (x <= m_xmin) return 0;
  if (x >= m_xmax) return m_nbins;

  if (m_mode == axis_mode::linear) {
    return static_cast<size_t>((x - m_xmin) / m_xstep);
  }
  // Assumes axis_mode::log
  return static_cast<size_t>(log(x / m_xmin) / m_xstep);
}

//////////////////////////////////////////////////////////////////////////////
// One-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
OneDim_Table::OneDim_Table(const axis& xbins) : m_x(xbins)
{
  if (m_x.m_nbins > 0) m_values.resize(m_x.m_nbins + 1, 0.0);
}

void OneDim_Table::Fill(size_t xbin, double value)
{
  if (xbin <= m_x.m_nbins) m_values[xbin] = value;
}

void OneDim_Table::Rescale(double factor)
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i) m_values[i] *= factor;
}

double OneDim_Table::operator()(double x) const
{
  if (x < m_x.m_xmin || x > m_x.m_xmax) return 0.0;

  size_t bin = m_x.bin(x);

  double x1 = m_x.x(bin);
  double dx = (m_x.m_mode == axis_mode::linear) ? m_x.m_xstep : m_x.x(bin + 1) - x1;

  if (dx < 1.e-12 * std::fabs(x1)) return m_values[bin];

  double w = (x - x1) / dx; // Normalized weight
  return m_values[bin] * (1.0 - w) + m_values[bin + 1] * w;
}

double OneDim_Table::Inverse(double value) const
{
  // This method only works for suitably normalised, monotonic "cumulative" tables.
  if (m_values.empty()) return 0.0;

  // Use binary search (std::lower_bound) for efficiency
  auto it = std::lower_bound(m_values.begin(), m_values.end(), value);

  if (it == m_values.begin()) return m_x.m_xmin;
  if (it == m_values.end()) return m_x.m_xmax;

  size_t bin = std::distance(m_values.begin(), it) - 1;

  double val1 = m_values[bin];
  double val2 = m_values[bin + 1];

  if (std::fabs(val2 - val1) < 1.e-12 * std::fabs(val1)) return m_x.x(bin);

  return (m_x.x(bin) * (val2 - value) + m_x.x(bin + 1) * (value - val1)) / (val2 - val1);
}

double OneDim_Table::Integral() const
{
  double integral = 0;
  for (size_t i = 0; i < m_x.m_nbins; ++i) {
    integral += (m_x.x(i + 1) - m_x.x(i)) * (m_values[i] + m_values[i + 1]) / 2.0;
  }
  return integral;
}

std::unique_ptr<OneDim_Table> OneDim_Table::Cumulative(double expo, double& integral) const
{
  auto cumulative = std::make_unique<OneDim_Table>(m_x);
  integral = 0.0;
  cumulative->Fill(0, 0.0);
  for (size_t i = 1; i <= m_x.m_nbins; ++i) {
    double estimate = 0;
    double x_mid = (m_x.x(i) + m_x.x(i - 1)) / 2.0;
    double dx = m_x.x(i) - m_x.x(i - 1);

    if (m_x.m_mode == axis_mode::linear)
      estimate = dx * (m_values[i] + m_values[i - 1]) / 2.0;
    else // axis_mode::log
      estimate = dx * std::sqrt(m_values[i] * m_values[i - 1]);

    integral += estimate * std::pow(x_mid, expo);
    cumulative->Fill(i, integral);
  }
  return cumulative;
}

void OneDim_Table::OutputToCSV(std::ofstream& outfile) const
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i) {
    outfile << m_x.x(i) << "," << m_values[i] << std::endl;
  }
}

//////////////////////////////////////////////////////////////////////////////
// Two-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
TwoDim_Table::TwoDim_Table(const axis& xbins, const axis& ybins) : m_x(xbins), m_y(ybins)
{
  if (m_y.m_nbins == 0) THROW(fatal_error, "Zero bins in y direction not supported.");
  m_values.resize((m_x.m_nbins + 1) * (m_y.m_nbins + 1), 0.0);
}

void TwoDim_Table::Fill(size_t xbin, size_t ybin, double value)
{
  if (xbin <= m_x.m_nbins && ybin <= m_y.m_nbins) m_values[Index(xbin, ybin)] = value;
}

double TwoDim_Table::operator()(double x, double y) const
{
  if (x < m_x.m_xmin || x > m_x.m_xmax || y < m_y.m_xmin || y > m_y.m_xmax) return 0.0;

  size_t xbin = m_x.bin(x);
  size_t ybin = m_y.bin(y);
  if (xbin >= m_x.m_nbins) xbin = m_x.m_nbins - 1;
  if (ybin >= m_y.m_nbins) ybin = m_y.m_nbins - 1;

  // Get weights for interpolation
  double x1 = m_x.x(xbin);
  double dx = (m_x.m_mode == axis_mode::linear) ? m_x.m_xstep : m_x.x(xbin + 1) - x1;
  double wx = (std::fabs(dx) > 0) ? (x - x1) / dx : 0.0;

  double y1 = m_y.x(ybin);
  double dy = (m_y.m_mode == axis_mode::linear) ? m_y.m_xstep : m_y.x(ybin + 1) - y1;
  double wy = (std::fabs(dy) > 0) ? (y - y1) / dy : 0.0;

  // Bilinear interpolation
  const double v00 = Value(xbin, ybin);
  const double v10 = Value(xbin + 1, ybin);
  const double v01 = Value(xbin, ybin + 1);
  const double v11 = Value(xbin + 1, ybin + 1);

  double v0 = v00 * (1.0 - wx) + v10 * wx;
  double v1 = v01 * (1.0 - wx) + v11 * wx;

  return v0 * (1.0 - wy) + v1 * wy;
}

double TwoDim_Table::Integral(size_t naxis, size_t bin) const
{
  double integral = 0.;
  if (naxis == 0) { // Integrate along X-axis for a fixed Y-bin
    if (bin > m_y.m_nbins) return 0.0;
    for (size_t i = 0; i < m_x.m_nbins; ++i)
      integral += (m_x.x(i + 1) - m_x.x(i)) * (Value(i + 1, bin) + Value(i, bin)) / 2.0;
  } else { // Integrate along Y-axis for a fixed X-bin
    if (bin > m_x.m_nbins) return 0.0;
    for (size_t i = 0; i < m_y.m_nbins; ++i)
      integral += (m_y.x(i + 1) - m_y.x(i)) * (Value(bin, i + 1) + Value(bin, i)) / 2.0;
  }
  return integral;
}

std::unique_ptr<OneDim_Table> TwoDim_Table::Cumulative(double expo, size_t naxis,
                                                       size_t bin, double& integral) const
{
  std::unique_ptr<OneDim_Table> cumulative;
  integral = 0.0;

  if (naxis == 0) { // Cumulative along X-axis
    if (bin > m_y.m_nbins) return std::make_unique<OneDim_Table>(m_x);
    cumulative = std::make_unique<OneDim_Table>(m_x);
    cumulative->Fill(0, 0.0);
    for (size_t i = 1; i <= m_x.m_nbins; ++i) {
      double estimate = (m_x.x(i) - m_x.x(i - 1)) * (Value(i, bin) + Value(i - 1, bin)) / 2.0;
      estimate *= std::pow((m_x.x(i) + m_x.x(i - 1)) / 2.0, expo);
      integral += estimate;
      cumulative->Fill(i, integral);
    }
  } else { // Cumulative along Y-axis
    if (bin > m_x.m_nbins) return std::make_unique<OneDim_Table>(m_y);
    cumulative = std::make_unique<OneDim_Table>(m_y);
    cumulative->Fill(0, 0.0);
    for (size_t i = 1; i <= m_y.m_nbins; ++i) {
      double estimate = (m_y.x(i) - m_y.x(i - 1)) * (Value(bin, i) + Value(bin, i - 1)) / 2.0;
      estimate *= std::pow((m_y.x(i) + m_y.x(i - 1)) / 2.0, expo);
      integral += estimate;
      cumulative->Fill(i, integral);
    }
  }
  return cumulative;
}

void TwoDim_Table::OutputToCSV(std::ofstream& outfile) const
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i)
    for (size_t j = 0; j <= m_y.m_nbins; ++j)
      outfile << m_x.x(i) << "," << m_y.x(j) << "," << Value(i, j) << std::endl;
}

//////////////////////////////////////////////////////////////////////////////
// Three-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
ThreeDim_Table::ThreeDim_Table(const axis& xbins, const axis& ybins, const axis& zbins)
    : m_x(xbins), m_y(ybins), m_z(zbins)
{
  if (m_z.m_nbins == 0) THROW(fatal_error, "Zero bins in z direction not supported.");
  m_values.resize((m_x.m_nbins + 1) * (m_y.m_nbins + 1) * (m_z.m_nbins + 1), 0.0);
}

void ThreeDim_Table::Fill(size_t xbin, size_t ybin, size_t zbin, double value)
{
  if (xbin <= m_x.m_nbins && ybin <= m_y.m_nbins && zbin <= m_z.m_nbins)
    m_values[Index(xbin, ybin, zbin)] = value;
}

inline double lerp(double v0, double v1, double w) {
  return v0 * (1.0 - w) + v1 * w;
}

double ThreeDim_Table::operator()(double x, double y, double z) const
{
  if (x < m_x.m_xmin || x > m_x.m_xmax || y < m_y.m_xmin || y > m_y.m_xmax ||
      z < m_z.m_xmin || z > m_z.m_xmax)
    return 0.;

  size_t xbin = m_x.bin(x);
  size_t ybin = m_y.bin(y);
  size_t zbin = m_z.bin(z);
  if (xbin >= m_x.m_nbins) xbin = m_x.m_nbins - 1;
  if (ybin >= m_y.m_nbins) ybin = m_y.m_nbins - 1;
  if (zbin >= m_z.m_nbins) zbin = m_z.m_nbins - 1;

  // Get weights for interpolation
  double x1 = m_x.x(xbin);
  double dx = (m_x.m_mode == axis_mode::linear) ? m_x.m_xstep : m_x.x(xbin + 1) - x1;
  double wx = (std::fabs(dx) > 0) ? (x - x1) / dx : 0.0;

  double y1 = m_y.x(ybin);
  double dy = (m_y.m_mode == axis_mode::linear) ? m_y.m_xstep : m_y.x(ybin + 1) - y1;
  double wy = (std::fabs(dy) > 0) ? (y - y1) / dy : 0.0;

  double z1 = m_z.x(zbin);
  double dz = (m_z.m_mode == axis_mode::linear) ? m_z.m_xstep : m_z.x(zbin + 1) - z1;
  double wz = (std::fabs(dz) > 0) ? (z - z1) / dz : 0.0;

  // Get 8 corner values
  const double v000 = Value(xbin, ybin, zbin);
  const double v100 = Value(xbin + 1, ybin, zbin);
  const double v010 = Value(xbin, ybin + 1, zbin);
  const double v110 = Value(xbin + 1, ybin + 1, zbin);
  const double v001 = Value(xbin, ybin, zbin + 1);
  const double v101 = Value(xbin + 1, ybin, zbin + 1);
  const double v011 = Value(xbin, ybin + 1, zbin + 1);
  const double v111 = Value(xbin + 1, ybin + 1, zbin + 1);

  // Trilinear interpolation
  double c00 = lerp(v000, v100, wx);
  double c10 = lerp(v010, v110, wx);
  double c01 = lerp(v001, v101, wx);
  double c11 = lerp(v011, v111, wx);

  double c0 = lerp(c00, c10, wy);
  double c1 = lerp(c01, c11, wy);

  return lerp(c0, c1, wz);
}

void ThreeDim_Table::OutputToCSV(std::ofstream& outfile) const
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i)
    for (size_t j = 0; j <= m_y.m_nbins; ++j)
      for (size_t k = 0; k <= m_z.m_nbins; ++k)
        outfile << m_x.x(i) << "," << m_y.x(j) << "," << m_z.x(k) << ","
                << Value(i,j,k) << std::endl;
}
