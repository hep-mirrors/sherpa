#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Org/Message.H"

#include <catch2/catch_all.hpp>
#include <cmath>
#include <sstream>

using namespace ATOOLS;

namespace {
  // THROW (the axis ctor's input checks) touches the message singleton, which
  // is null in a bare test binary. Create it once.
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

TEST_CASE("TwoDim_Table interpolates in its 2nd axis",
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

  SECTION("single-bin x-axis (fixed energy) is flat in x") {
    const double S = 13000. * 13000.;
    axis ax(1, S, (1. + 1.e-6) * S, axis_mode::linear);
    axis ay(4, 0., 1., axis_mode::linear);
    TwoDim_Table t(ax, ay);
    for (size_t j = 0; j <= ay.m_nbins; ++j) t.Fill(0, j, ay.x(j)); // value = y
    CHECK(t(S, 0.3) == Catch::Approx(0.3));
    CHECK(t(S, 0.75) == Catch::Approx(0.75));
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
