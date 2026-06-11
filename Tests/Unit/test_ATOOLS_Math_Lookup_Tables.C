#include "ATOOLS/Math/Lookup_Tables.H"

#include <catch2/catch_all.hpp>
#include <cmath>
#include <sstream>

using namespace ATOOLS;

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
