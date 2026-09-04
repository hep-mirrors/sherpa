#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm> // For std::lower_bound, std::min
#include <cmath>     // For std::fabs, std::sqrt, std::exp, etc.
#include <fstream>
#include <iostream>

using namespace ATOOLS;

namespace {
  // Text helpers for table (de)serialisation. A text format is portable across
  // builds/endianness and round-trips through ATOOLS::My_File, which is how
  // Sherpa shares result files across MPI ranks (rank 0 reads the file and
  // broadcasts its content). Doubles are written with 17 significant digits
  // (std::numeric_limits<double>::max_digits10), which round-trips an IEEE
  // double exactly, so a reloaded table reproduces a freshly computed one.

  // Sanity cap on the bin count read from a (regenerable) cache file, so a
  // corrupt count cannot trigger a runaway allocation; any real table is far
  // smaller, and a violation is treated as a corrupt file (i.e. recompute).
  constexpr size_t s_max_bins = 100000000; // 1e8

  // The per-axis cap above still permits products like nx*ny ~ 1e16 for a 2D
  // table, so bound the total element count too (checked without overflowing
  // size_t). Caps the allocation at ~0.8 GB; any real table is far smaller.
  constexpr size_t s_max_elements = 100000000; // 1e8 doubles

  void WriteAxis(std::ostream& os, const axis& a)
  {
    os << a.m_nbins << ' ' << a.m_xmin << ' ' << a.m_xmax << ' '
       << static_cast<int>(a.m_mode) << '\n';
  }
  axis ReadAxis(std::istream& is)
  {
    size_t nbins = 0;
    double xmin = 0., xmax = 0.;
    int mode = 0;
    is >> nbins >> xmin >> xmax >> mode;
    if (!is || nbins < 1 || nbins > s_max_bins ||
        (mode != axis_mode::linear && mode != axis_mode::log))
      THROW(fatal_error, "Corrupt axis while reading cache.");
    return axis(nbins, xmin, xmax, static_cast<axis_mode::code>(mode));
  }
  void WriteValues(std::ostream& os, const std::vector<double>& v)
  {
    for (double x : v)
      os << x << '\n';
  }
  // Fills an already-sized vector (its length is fixed by the axes), so no
  // length is read from the file.
  void ReadValues(std::istream& is, std::vector<double>& v)
  {
    for (double& x : v)
      is >> x;
  }
} // namespace

axis::axis(size_t nbins, double xmin, double xmax, axis_mode::code mode)
    : m_nbins(nbins), m_xmin(xmin), m_xmax(xmax), m_mode(mode), m_xstep(0.0)
{
  if (m_nbins < 1) THROW(fatal_error, "Cannot create table with no bins.");
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

double axis::weight(size_t bin, double x) const
{
  // A single-bin axis is flat (only index 0 is ever filled) and a degenerate
  // one has no width; in both cases the caller must fall back on the stored
  // value, i.e. weight zero, rather than scale it by a spurious factor.
  if (m_nbins == 1 || m_xstep <= 0.) return 0.;
  if (bin >= m_nbins) bin = m_nbins - 1;
  if (m_mode == axis_mode::linear) return (x - m_xmin) / m_xstep - double(bin);
  // Log axis: the nodes are spaced geometrically, x_i = x_min exp(i * step),
  // so the interpolation fraction has to be taken in log(x). Weighting
  // linearly in x between geometrically spaced nodes systematically
  // over-weights the upper node -- an O(step^2) bias on a smooth function,
  // which is at the few-per-mille level for the ~20 nodes per decade the EPA
  // flux tables use. This mirrors bin() above, so the returned weight is in
  // [0, 1] whenever bin == this->bin(x) (and exactly 1 at x = x_max, where
  // the callers clamp the bin index down by one).
  return std::log(x / m_xmin) / m_xstep - double(bin);
}

double axis::coordinate(size_t bin, double weight) const
{
  // Keep the same flat-axis convention as weight(). Also preserve x()'s
  // upper-edge behaviour for inverse tables whose requested value lies beyond
  // a particular row's range.
  if (m_nbins == 1 || m_xstep <= 0. || bin >= m_nbins) return x(bin);
  if (m_mode == axis_mode::linear) return x(bin) + weight * m_xstep;
  return x(bin) * std::exp(weight * m_xstep);
}

//////////////////////////////////////////////////////////////////////////////
// One-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
OneDim_Table::OneDim_Table(const axis& xbins): m_x(xbins)
{
  if (m_x.m_nbins > 0) m_values.resize(m_x.m_nbins + 1, 0.0);
}

void OneDim_Table::Fill(size_t xbin, double value)
{
  if (xbin <= m_x.m_nbins) m_values[xbin] = value;
}

void OneDim_Table::Rescale(double factor)
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i)
    m_values[i] *= factor;
}

double OneDim_Table::operator()(double x) const
{
  if (x < m_x.m_xmin || x > m_x.m_xmax) return 0.0;
  // A single-bin table is flat: only m_values[0] is filled, so interpolating
  // would scale it by a spurious weight (x may sit far from the lone node).
  // Return the one stored value.
  if (m_x.m_nbins == 1) return m_values[0];

  size_t bin = m_x.bin(x);
  // bin() returns m_nbins at the upper edge, where m_values[bin + 1] would be
  // one past the end; clamp as the 2D/3D lookups already do.
  if (bin >= m_x.m_nbins) bin = m_x.m_nbins - 1;

  const double w = m_x.weight(bin, x);
  return m_values[bin] * (1.0 - w) + m_values[bin + 1] * w;
}

double OneDim_Table::Inverse(double value) const
{
  // This method only works for suitably normalised, monotonic "cumulative"
  // tables.
  if (m_values.empty()) return 0.0;

  // Use binary search (std::lower_bound) for efficiency
  auto it = std::lower_bound(m_values.begin(), m_values.end(), value);

  if (it == m_values.begin()) return m_x.m_xmin;
  if (it == m_values.end()) return m_x.m_xmax;

  size_t bin = std::distance(m_values.begin(), it) - 1;

  double val1 = m_values[bin];
  double val2 = m_values[bin + 1];

  if (std::fabs(val2 - val1) < 1.e-12 * std::fabs(val1)) return m_x.x(bin);

  const double weight = (value - val1) / (val2 - val1);
  return m_x.coordinate(bin, weight);
}

double OneDim_Table::Integral() const
{
  double integral = 0;
  for (size_t i = 0; i < m_x.m_nbins; ++i) {
    integral +=
        (m_x.x(i + 1) - m_x.x(i)) * (m_values[i] + m_values[i + 1]) / 2.0;
  }
  return integral;
}

std::unique_ptr<OneDim_Table> OneDim_Table::Cumulative(double expo,
                                                       double& integral) const
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

void OneDim_Table::Write(std::ostream& os) const
{
  os.precision(17);
  WriteAxis(os, m_x);
  WriteValues(os, m_values);
}

std::unique_ptr<OneDim_Table> OneDim_Table::Read(std::istream& is)
{
  const axis x = ReadAxis(is);
  auto table = std::make_unique<OneDim_Table>(x);
  ReadValues(is, table->m_values);
  if (!is)
    THROW(fatal_error, "Corrupt OneDim_Table while reading cache.");
  return table;
}

//////////////////////////////////////////////////////////////////////////////
// Two-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
TwoDim_Table::TwoDim_Table(const axis& xbins, const axis& ybins)
    : m_x(xbins), m_y(ybins)
{
  if (m_y.m_nbins == 0)
    THROW(fatal_error, "Zero bins in y direction not supported.");
  m_values.resize((m_x.m_nbins + 1) * (m_y.m_nbins + 1), 0.0);
}

void TwoDim_Table::Fill(size_t xbin, size_t ybin, double value)
{
  if (xbin <= m_x.m_nbins && ybin <= m_y.m_nbins)
    m_values[Index(xbin, ybin)] = value;
}

double TwoDim_Table::operator()(double x, double y) const
{
  if (x < m_x.m_xmin || x > m_x.m_xmax || y < m_y.m_xmin || y > m_y.m_xmax)
    return 0.0;

  size_t xbin = m_x.bin(x);
  size_t ybin = m_y.bin(y);
  if (xbin >= m_x.m_nbins) xbin = m_x.m_nbins - 1;
  if (ybin >= m_y.m_nbins) ybin = m_y.m_nbins - 1;

  // Interpolation weights in each axis' own metric (see axis::weight), which
  // also handles the flat single-bin and degenerate cases.
  const double wx = m_x.weight(xbin, x);
  const double wy = m_y.weight(ybin, y);

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
      integral +=
          (m_x.x(i + 1) - m_x.x(i)) * (Value(i + 1, bin) + Value(i, bin)) / 2.0;
  } else { // Integrate along Y-axis for a fixed X-bin
    if (bin > m_x.m_nbins) return 0.0;
    for (size_t i = 0; i < m_y.m_nbins; ++i)
      integral +=
          (m_y.x(i + 1) - m_y.x(i)) * (Value(bin, i + 1) + Value(bin, i)) / 2.0;
  }
  return integral;
}

std::unique_ptr<OneDim_Table> TwoDim_Table::Cumulative(double expo,
                                                       size_t naxis, size_t bin,
                                                       double& integral) const
{
  std::unique_ptr<OneDim_Table> cumulative;
  integral = 0.0;

  if (naxis == 0) { // Cumulative along X-axis
    if (bin > m_y.m_nbins) return std::make_unique<OneDim_Table>(m_x);
    cumulative = std::make_unique<OneDim_Table>(m_x);
    cumulative->Fill(0, 0.0);
    for (size_t i = 1; i <= m_x.m_nbins; ++i) {
      double estimate =
          (m_x.x(i) - m_x.x(i - 1)) * (Value(i, bin) + Value(i - 1, bin)) / 2.0;
      estimate *= std::pow((m_x.x(i) + m_x.x(i - 1)) / 2.0, expo);
      integral += estimate;
      cumulative->Fill(i, integral);
    }
  } else { // Cumulative along Y-axis
    if (bin > m_x.m_nbins) return std::make_unique<OneDim_Table>(m_y);
    cumulative = std::make_unique<OneDim_Table>(m_y);
    cumulative->Fill(0, 0.0);
    for (size_t i = 1; i <= m_y.m_nbins; ++i) {
      double estimate =
          (m_y.x(i) - m_y.x(i - 1)) * (Value(bin, i) + Value(bin, i - 1)) / 2.0;
      estimate *= std::pow((m_y.x(i) + m_y.x(i - 1)) / 2.0, expo);
      integral += estimate;
      cumulative->Fill(i, integral);
    }
  }
  return cumulative;
}

std::unique_ptr<TwoDim_Table> TwoDim_Table::Invert(size_t axislabel,
                                                   size_t nbins) const
{
  if (axislabel != 1)
    THROW(fatal_error, "Two-dim table inversion only available for 2nd axis.");
  // Determine the monotonicity direction and the global value range over all
  // x rows. The last y grid point is index m_y.m_nbins (intervals convention).
  const size_t ylast = m_y.m_nbins;
  bool   dir    = Value(0, 0) > Value(0, ylast);
  double minval = dir ? Value(0, ylast) : Value(0, 0);
  double maxval = dir ? Value(0, 0) : Value(0, ylast);
  for (size_t i = 1; i <= m_x.m_nbins; ++i) {
    if (dir) {
      if (Value(i, ylast) < minval) minval = Value(i, ylast);
      if (Value(i, 0)     > maxval) maxval = Value(i, 0);
    } else {
      if (Value(i, ylast) > maxval) maxval = Value(i, ylast);
      if (Value(i, 0)     < minval) minval = Value(i, 0);
    }
  }
  axis vaxis(nbins, minval, maxval, axis_mode::linear);
  auto table = std::make_unique<TwoDim_Table>(m_x, vaxis);
  for (size_t i = 0; i <= m_x.m_nbins; ++i) {
    for (size_t j = 0; j <= nbins; ++j) {
      double value = vaxis.x(j), val1 = 0., val2 = 0.;
      size_t ybin = 0;
      for (ybin = 0; ybin < m_y.m_nbins; ++ybin) {
        val1 = Value(i, ybin);
        val2 = Value(i, ybin + 1);
        if (dir) {
          if (val1 >= value && val2 <= value) break;
        } else {
          if (val1 <= value && val2 >= value) break;
        }
      }
      // The same expression gives a weight in [0,1] for increasing and
      // decreasing rows. Convert it back through the axis metric so this is the
      // actual inverse of operator() on both linear and logarithmic axes.
      const double weight = val1 != val2 ? (value - val1) / (val2 - val1) : 0.;
      table->Fill(i, j, m_y.coordinate(ybin, weight));
    }
  }
  return table;
}

void TwoDim_Table::OutputToCSV(std::ofstream& outfile) const
{
  for (size_t i = 0; i <= m_x.m_nbins; ++i)
    for (size_t j = 0; j <= m_y.m_nbins; ++j)
      outfile << m_x.x(i) << "," << m_y.x(j) << "," << Value(i, j) << std::endl;
}

void TwoDim_Table::Write(std::ostream& os) const
{
  os.precision(17);
  WriteAxis(os, m_x);
  WriteAxis(os, m_y);
  WriteValues(os, m_values);
}

std::unique_ptr<TwoDim_Table> TwoDim_Table::Read(std::istream& is)
{
  const axis x = ReadAxis(is);
  const axis y = ReadAxis(is);
  // ReadAxis guarantees nbins >= 1, so ny >= 2 and the division is safe; this
  // bounds the product nx*ny without overflowing size_t.
  const size_t nx = x.m_nbins + 1, ny = y.m_nbins + 1;
  if (nx > s_max_elements / ny)
    THROW(fatal_error, "Corrupt TwoDim_Table size while reading cache.");
  auto table = std::make_unique<TwoDim_Table>(x, y);
  ReadValues(is, table->m_values);
  if (!is)
    THROW(fatal_error, "Corrupt TwoDim_Table while reading cache.");
  return table;
}

//////////////////////////////////////////////////////////////////////////////
// Three-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
ThreeDim_Table::ThreeDim_Table(const axis& xbins, const axis& ybins,
                               const axis& zbins)
    : m_x(xbins), m_y(ybins), m_z(zbins)
{
  if (m_z.m_nbins == 0)
    THROW(fatal_error, "Zero bins in z direction not supported.");
  m_values.resize((m_x.m_nbins + 1) * (m_y.m_nbins + 1) * (m_z.m_nbins + 1),
                  0.0);
}

void ThreeDim_Table::Fill(size_t xbin, size_t ybin, size_t zbin, double value)
{
  if (xbin <= m_x.m_nbins && ybin <= m_y.m_nbins && zbin <= m_z.m_nbins)
    m_values[Index(xbin, ybin, zbin)] = value;
}

inline double lerp(double v0, double v1, double w)
{
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

  // Interpolation weights in each axis' own metric (see axis::weight), which
  // also handles the flat single-bin and degenerate cases.
  const double wx = m_x.weight(xbin, x);
  const double wy = m_y.weight(ybin, y);
  const double wz = m_z.weight(zbin, z);

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
                << Value(i, j, k) << std::endl;
}
