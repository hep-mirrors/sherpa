#include <catch2/catch_all.hpp>

#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Org/Message.H"

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
