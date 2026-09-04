#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Org/Message.H"

#include <catch2/catch_all.hpp>
#include <cmath>
#include <sstream>

using namespace ATOOLS;

namespace {
  // THROW (the axis ctor's input checks and Invert's axislabel guard) touches
  // the message singleton, which is null in a bare test binary. Create it once.
  void EnsureMsg() { if (!msg) msg = new Message(); }
}

TEST_CASE("axis uses the intervals convention", "[ATOOLS::axis]") {
  EnsureMsg();

  SECTION("linear: nbins is the interval count, x(nbins) is the upper edge") {
    axis a(4, 0., 4., axis_mode::linear);
    CHECK(a.m_nbins == 4);
    CHECK(a.m_xstep == Catch::Approx(1.0));   // (4-0)/4, NOT /(nbins-1)
    CHECK(a.x(0) == Catch::Approx(0.0));
    CHECK(a.x(2) == Catch::Approx(2.0));
    CHECK(a.x(4) == Catch::Approx(4.0));      // nbins+1 grid points: 0..4
    CHECK(a.bin(0.0) == 0);
    CHECK(a.bin(2.5) == 2);
    CHECK(a.bin(4.0) == 4);                   // upper edge snaps to nbins
  }

  SECTION("log: geometric spacing") {
    axis a(2, 1., 100., axis_mode::log);
    CHECK(a.x(0) == Catch::Approx(1.0));
    CHECK(a.x(1) == Catch::Approx(10.0));
    CHECK(a.x(2) == Catch::Approx(100.0));
  }
}

TEST_CASE("OneDim_Table interpolates and inverts a cumulative",
          "[ATOOLS::OneDim_Table]") {
  EnsureMsg();

  SECTION("linear interpolation of a linear function is exact") {
    axis a(4, 0., 4., axis_mode::linear);
    OneDim_Table t(a);
    for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, a.x(i)); // v(x) = x
    CHECK(t(1.5) == Catch::Approx(1.5));
    CHECK(t(3.0) == Catch::Approx(3.0));
    CHECK(t(0.0) == Catch::Approx(0.0));
  }

  SECTION("Cumulative + Inverse round-trip on a flat density") {
    axis a(4, 0., 4., axis_mode::linear);
    OneDim_Table t(a);
    for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, 1.0); // density == 1
    double integral = 0.;
    auto cum = t.Cumulative(0.0, integral);
    CHECK(integral == Catch::Approx(4.0));     // area under unit density over [0,4]
    CHECK((*cum)(2.0) == Catch::Approx(2.0));  // cumulative(x) == x
    CHECK(cum->Inverse(2.0) == Catch::Approx(2.0));
    CHECK(cum->Inverse(3.5) == Catch::Approx(3.5));
  }

  SECTION("Inverse follows the metric of a logarithmic axis") {
    axis a(4, 1., 100., axis_mode::log);
    OneDim_Table t(a);
    for (size_t i = 0; i <= a.m_nbins; ++i)
      t.Fill(i, std::log(a.x(i)));
    for (double x : {1.7, 4.2, 13., 57.})
      CHECK(t.Inverse(t(x)) == Catch::Approx(x).epsilon(1.e-12));
  }

  SECTION("Rescale scales every entry") {
    axis a(2, 0., 2., axis_mode::linear);
    OneDim_Table t(a);
    for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, 2.0);
    t.Rescale(3.0);
    CHECK(t(1.0) == Catch::Approx(6.0));
  }

  SECTION("single-bin (fixed-energy) table is flat") {
    // Mirrors the fixed-cms s-axis: axis(1, S, (1+1e-6)S). With a large S a
    // naive interpolation weight blows up; the table must return the one value.
    const double S = 13000. * 13000.;
    axis a(1, S, (1. + 1.e-6) * S, axis_mode::linear);
    OneDim_Table t(a);
    t.Fill(0, 42.0);
    CHECK(t(S) == Catch::Approx(42.0));
    CHECK(t((1. + 0.5e-6) * S) == Catch::Approx(42.0));
  }
}

TEST_CASE("TwoDim_Table interpolates and inverts in its 2nd axis",
          "[ATOOLS::TwoDim_Table]") {
  EnsureMsg();

  SECTION("bilinear interpolation of x+y is exact") {
    axis ax(2, 0., 2., axis_mode::linear), ay(2, 0., 2., axis_mode::linear);
    TwoDim_Table t(ax, ay);
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j) t.Fill(i, j, ax.x(i) + ay.x(j));
    CHECK(t(0.5, 0.5) == Catch::Approx(1.0));
    CHECK(t(1.5, 0.5) == Catch::Approx(2.0));
  }

  SECTION("Invert round-trips an increasing CDF (y), v -> y") {
    axis ax(2, 0., 1., axis_mode::linear), ay(4, 0., 1., axis_mode::linear);
    TwoDim_Table cdf(ax, ay);
    // T(x,y) = y : increasing in y, independent of x -> inverse is I(x,v) = v.
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j) cdf.Fill(i, j, ay.x(j));
    auto inv = cdf.Invert(1, 20);
    CHECK((*inv)(0.5, 0.0) == Catch::Approx(0.0).margin(1e-9));
    CHECK((*inv)(0.5, 0.3) == Catch::Approx(0.3));
    CHECK((*inv)(0.5, 1.0) == Catch::Approx(1.0));
  }

  SECTION("Invert round-trips a decreasing CDF (y), v -> 1-y") {
    axis ax(2, 0., 1., axis_mode::linear), ay(4, 0., 1., axis_mode::linear);
    TwoDim_Table cdf(ax, ay);
    // T(x,y) = 1-y : decreasing in y -> inverse is I(x,v) = 1-v.
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j) cdf.Fill(i, j, 1.0 - ay.x(j));
    auto inv = cdf.Invert(1, 20);
    CHECK((*inv)(0.5, 0.3) == Catch::Approx(0.7));
    CHECK((*inv)(0.5, 0.8) == Catch::Approx(0.2));
  }

  SECTION("Invert follows the metric of a logarithmic y-axis") {
    axis ax(2, 0., 1., axis_mode::linear), ay(4, 1., 100., axis_mode::log);
    TwoDim_Table cdf(ax, ay);
    // T(x,y) = log(y)/log(100), so its inverse is y = 100^v. Use eight
    // inverse-table intervals to probe both nodes and midpoints of the input
    // y-axis while staying on output-table nodes.
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j)
        cdf.Fill(i, j, std::log(ay.x(j)) / std::log(100.));
    auto inv = cdf.Invert(1, 8);
    for (size_t j = 0; j <= 8; ++j) {
      const double v = double(j) / 8.;
      CHECK((*inv)(0.5, v) == Catch::Approx(std::pow(100., v)).epsilon(1.e-12));
    }
  }

  SECTION("Invert rejects inversion of the 1st axis") {
    axis ax(2, 0., 1., axis_mode::linear), ay(2, 0., 1., axis_mode::linear);
    TwoDim_Table t(ax, ay);
    REQUIRE_THROWS(t.Invert(0, 10));
  }

  SECTION("single-bin x-axis (fixed energy) is flat in x") {
    const double S = 13000. * 13000.;
    axis ax(1, S, (1. + 1.e-6) * S, axis_mode::linear);
    axis ay(4, 0., 1., axis_mode::linear);
    TwoDim_Table t(ax, ay);
    for (size_t j = 0; j <= ay.m_nbins; ++j) t.Fill(0, j, ay.x(j)); // value = y
    CHECK(t(S, 0.3) == Catch::Approx(0.3));
    CHECK(t(S, 0.75) == Catch::Approx(0.75));
  }

  SECTION("Invert spans the full range when the CDF reaches 1 at the top node") {
    // Regression guard for Interaction_Probability::FillIntegrated: the CDF must
    // be filled across the WHOLE node grid (0..m_nbins) so it reaches 1.0 at the
    // top node that Invert reads as the maximum. A consumer that fills only
    // < m_nbins leaves that node at 0, collapsing the inverted value-axis to
    // [0,0] and making SelectB return ~0 for every draw.
    axis ax(2, 0., 1., axis_mode::linear), ay(8, 0., 4., axis_mode::linear);
    TwoDim_Table cdf(ax, ay);
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j)
        cdf.Fill(i, j, ay.x(j) / ay.m_xmax); // 0 at b=0, exactly 1 at b=bmax
    auto inv = cdf.Invert(1, 50);
    CHECK((*inv)(0.5, 0.0) == Catch::Approx(0.0).margin(1e-9));
    CHECK((*inv)(0.5, 0.5) == Catch::Approx(2.0));
    CHECK((*inv)(0.5, 1.0) == Catch::Approx(4.0)); // full b-range, not collapsed
  }
}

TEST_CASE("OneDim_Table survives a binary round-trip",
          "[ATOOLS::Lookup_Tables]")
{
  axis a(50, 0.1, 10., axis_mode::log);
  OneDim_Table t(a);
  for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, std::sin(0.3 * i) + i);

  std::stringstream ss; // text format, round-trips exactly at 17 sig. digits
  t.Write(ss);
  auto t2 = OneDim_Table::Read(ss);
  REQUIRE(t2);

  SECTION("axis and stored values are reproduced exactly")
  {
    CHECK(t2->NBins() == a.m_nbins);
    CHECK(t2->GetAxis().m_xmin == a.m_xmin);
    CHECK(t2->GetAxis().m_xmax == a.m_xmax);
    CHECK(t2->GetAxis().m_mode == a.m_mode);
    for (size_t i = 0; i <= a.m_nbins; ++i) CHECK(t2->Value(i) == t.Value(i));
  }

  SECTION("interpolation agrees at intermediate points")
  {
    for (double x : {0.2, 1.0, 3.3, 7.7})
      CHECK((*t2)(x) == Catch::Approx(t(x)));
  }
}

TEST_CASE("TwoDim_Table survives a binary round-trip",
          "[ATOOLS::Lookup_Tables]")
{
  axis ax(20, 1.e-3, 1., axis_mode::log);
  axis ay(15, 0., 5., axis_mode::linear);
  TwoDim_Table t(ax, ay);
  for (size_t i = 0; i <= ax.m_nbins; ++i)
    for (size_t j = 0; j <= ay.m_nbins; ++j)
      t.Fill(i, j, 0.5 * i - j + 0.1 * i * j);

  std::stringstream ss; // text format, round-trips exactly at 17 sig. digits
  t.Write(ss);
  auto t2 = TwoDim_Table::Read(ss);
  REQUIRE(t2);

  SECTION("stored values are reproduced exactly")
  {
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j)
        CHECK(t2->Value(i, j) == t.Value(i, j));
  }

  SECTION("bilinear interpolation agrees at intermediate points")
  {
    for (double x : {1.e-2, 0.1, 0.5})
      for (double y : {0.5, 2.0, 4.0})
        CHECK((*t2)(x, y) == Catch::Approx(t(x, y)));
  }
}

////////////////////////////////////////////////////////////////////////////////
// Interpolation on a log axis happens in log(x).
//
// The nodes of a log axis are spaced geometrically, x_i = x_min exp(i*step), so
// the interpolation fraction has to be taken in log(x). Weighting linearly in x
// between geometrically spaced nodes systematically over-weights the upper node
// -- an O(step^2) bias on every off-node lookup, a few per mille at the ~20
// nodes per decade the EPA flux tables use.
//
// The sharp characterisation used below: with the log metric, a function that
// is *linear in log(x)* is reproduced exactly at arbitrary off-node points.
// Nothing above catches this -- the round-trip tests compare a table against
// itself, so a shared systematic bias cancels.
////////////////////////////////////////////////////////////////////////////////

TEST_CASE("log axes interpolate in log(x)", "[ATOOLS::axis][ATOOLS::OneDim_Table]")
{
  EnsureMsg();

  SECTION("axis::weight is the fraction in log(x), 1/2 at the geometric mean")
  {
    axis a(2, 1., 100., axis_mode::log); // nodes 1, 10, 100
    CHECK(a.weight(0, 1.0) == Catch::Approx(0.0));
    CHECK(a.weight(0, 10.0) == Catch::Approx(1.0));
    // sqrt(1*10) is the *geometric* midpoint of the first bin
    CHECK(a.weight(0, std::sqrt(10.0)) == Catch::Approx(0.5));
    CHECK(a.weight(1, std::sqrt(1000.0)) == Catch::Approx(0.5));
    CHECK(a.coordinate(0, 0.5) == Catch::Approx(std::sqrt(10.0)));
    CHECK(a.coordinate(1, 0.5) == Catch::Approx(std::sqrt(1000.0)));
  }

  SECTION("linear axes are unaffected: fraction is still the fraction in x")
  {
    axis a(4, 0., 4., axis_mode::linear);
    CHECK(a.weight(0, 0.0) == Catch::Approx(0.0));
    CHECK(a.weight(2, 2.5) == Catch::Approx(0.5));
    CHECK(a.weight(3, 4.0) == Catch::Approx(1.0));
    CHECK(a.coordinate(2, 0.5) == Catch::Approx(2.5));
  }

  SECTION("OneDim_Table: a function linear in log(x) is exact off-node")
  {
    axis a(20, 1.e-3, 1., axis_mode::log);
    OneDim_Table t(a);
    const auto f = [](double x) { return 2.5 + 0.75 * std::log(x); };
    for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, f(a.x(i)));
    for (double x : {1.3e-3, 7.7e-3, 4.2e-2, 0.31, 0.94})
      CHECK(t(x) == Catch::Approx(f(x)).epsilon(1.e-12));
  }

  SECTION("TwoDim_Table: bilinear in (log x, log y) is exact off-node")
  {
    axis ax(12, 1.e-2, 1.e2, axis_mode::log);
    axis ay(8, 1.e-1, 1.e1, axis_mode::log);
    TwoDim_Table t(ax, ay);
    const auto f = [](double x, double y) {
      return 1.5 + 0.5 * std::log(x) - 0.25 * std::log(y);
    };
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j)
        t.Fill(i, j, f(ax.x(i), ay.x(j)));
    for (double x : {3.3e-2, 1.7, 47.})
      for (double y : {0.17, 1.9, 7.3})
        CHECK(t(x, y) == Catch::Approx(f(x, y)).epsilon(1.e-12));
  }

  SECTION("ThreeDim_Table: trilinear in (log x, log y, log z) is exact off-node")
  {
    axis ax(6, 1.e-2, 1.e2, axis_mode::log);
    axis ay(6, 1.e-1, 1.e1, axis_mode::log);
    axis az(6, 1., 1.e3, axis_mode::log);
    ThreeDim_Table t(ax, ay, az);
    const auto f = [](double x, double y, double z) {
      return 0.5 * std::log(x) - 0.25 * std::log(y) + 0.125 * std::log(z);
    };
    for (size_t i = 0; i <= ax.m_nbins; ++i)
      for (size_t j = 0; j <= ay.m_nbins; ++j)
        for (size_t k = 0; k <= az.m_nbins; ++k)
          t.Fill(i, j, k, f(ax.x(i), ay.x(j), az.x(k)));
    for (double x : {3.3e-2, 47.})
      for (double y : {0.17, 7.3})
        for (double z : {6.1, 410.})
          CHECK(t(x, y, z) == Catch::Approx(f(x, y, z)).epsilon(1.e-12));
  }

  SECTION("upper edge of a log axis returns the last node, not past it")
  {
    // bin() returns m_nbins at x_max, where values[bin+1] would be one past
    // the end; the lookups clamp the bin index and the weight becomes 1.
    axis a(5, 1.e-2, 1., axis_mode::log);
    OneDim_Table t(a);
    for (size_t i = 0; i <= a.m_nbins; ++i) t.Fill(i, double(i));
    CHECK(t(1.0) == Catch::Approx(5.0));
    CHECK(a.weight(a.m_nbins - 1, 1.0) == Catch::Approx(1.0));
  }
}
