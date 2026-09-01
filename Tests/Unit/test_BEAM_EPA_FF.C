#include "BEAM/Spectra/EPA_FF.H"

#include <catch2/catch_all.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace BEAM;

namespace {
  // Independent restatement of the upper sampling limit used by
  // EPA_FF_Base::SampleB, in dimensionless units of R.
  double bmax_dimless(double x, double bmin, double bmax, double chi_cut,
                      double mass, double R)
  {
    return std::max(bmin, std::min(bmax, chi_cut / (x * mass * R)));
  }
}

TEST_CASE("EPA_FF_Base::SampleB caps the b-sampling at the 1/x flux support",
          "[BEAM::EPA_FF]")
{
  // representative proton-like parameters (dimensionless b in units of R)
  const double bmin = 0.3, bmax = 1.e3, chi_cut = 20., mass = 0.938, R = 3.5;

  SECTION("sampled b stays within [bmin, b_max(x)] for all ran")
  {
    for (double x : {1.e-4, 1.e-2, 0.1, 0.5, 0.9}) {
      const double bmx = bmax_dimless(x, bmin, bmax, chi_cut, mass, R);
      for (int i = 0; i <= 100; ++i) {
        double ran = std::min(i / 100.0, 1.0 - 1e-9);
        double b_phys = -1.;
        EPA_FF_Base::SampleB(x, ran, bmin, bmax, chi_cut, mass, R, b_phys);
        INFO("x=" << x << " ran=" << ran << " b_phys=" << b_phys
                  << " b_max_phys=" << bmx * R);
        CHECK(b_phys >= bmin * R * (1. - 1e-9));
        CHECK(b_phys <= bmx * R * (1. + 1e-9));
      }
    }
  }

  SECTION("endpoints: ran=0 -> bmin*R, ran->1 -> b_max(x)*R")
  {
    const double x = 0.5;
    const double bmx = bmax_dimless(x, bmin, bmax, chi_cut, mass, R);
    double b0 = -1., b1 = -1.;
    EPA_FF_Base::SampleB(x, 0.0, bmin, bmax, chi_cut, mass, R, b0);
    EPA_FF_Base::SampleB(x, 1.0 - 1e-12, bmin, bmax, chi_cut, mass, R, b1);
    CHECK_THAT(b0, Catch::Matchers::WithinRel(bmin * R, 1e-6));
    CHECK_THAT(b1, Catch::Matchers::WithinRel(bmx * R, 1e-6));
  }

  SECTION("support scales as 1/x where the cap is active")
  {
    // both x values give chi_cut/(x m R) inside (bmin, bmax) -> cap active
    const double x1 = 0.2, x2 = 0.4;
    const double b1 = bmax_dimless(x1, bmin, bmax, chi_cut, mass, R);
    const double b2 = bmax_dimless(x2, bmin, bmax, chi_cut, mass, R);
    CHECK_THAT(b1 / b2, Catch::Matchers::WithinRel(x2 / x1, 1e-9));
  }

  SECTION("clamps: small x -> bMax, large x -> bMin")
  {
    CHECK_THAT(bmax_dimless(1.e-8, bmin, bmax, chi_cut, mass, R),
               Catch::Matchers::WithinRel(bmax, 1e-12));
    CHECK_THAT(bmax_dimless(1.e3, bmin, bmax, chi_cut, mass, R),
               Catch::Matchers::WithinRel(bmin, 1e-12));
  }

  SECTION("importance-sampling identity: E_ran[weight] = (b_max(x)-bmin)*R")
  {
    // The proposal pdf integrates to 1, so the mean returned weight over a
    // uniform ran must equal the physical sampling range (b_max(x)-bmin)*R.
    for (double x : {1.e-3, 0.1, 0.7}) {
      const double bmx = bmax_dimless(x, bmin, bmax, chi_cut, mass, R);
      const int N = 200000;
      double sum = 0.;
      for (int i = 0; i < N; ++i) {
        const double ran = (i + 0.5) / N; // midpoint quadrature of E_ran
        double b_phys = -1.;
        sum +=
            EPA_FF_Base::SampleB(x, ran, bmin, bmax, chi_cut, mass, R, b_phys);
      }
      const double mean = sum / N;
      const double expected = (bmx - bmin) * R;
      INFO("x=" << x << " mean=" << mean << " expected=" << expected);
      CHECK_THAT(mean, Catch::Matchers::WithinRel(expected, 1e-3));
    }
  }
}

TEST_CASE("EPA_ff_type round-trips through operator<< / operator>>",
          "[BEAM::EPA_FF]")
{
  const std::vector<EPA_ff_type> all{
      EPA_ff_type::point,        EPA_ff_type::pointApprox,
      EPA_ff_type::proton,       EPA_ff_type::protonApprox,
      EPA_ff_type::ionApproxInt, EPA_ff_type::Gauss,
      EPA_ff_type::hcs,          EPA_ff_type::dipole,
      EPA_ff_type::dipoleApprox, EPA_ff_type::WoodSaxon,
      EPA_ff_type::WoodSaxonApprox, EPA_ff_type::ionApprox,
      EPA_ff_type::Test,         EPA_ff_type::Undefined};

  SECTION("every value survives a << then >> round-trip")
  {
    for (EPA_ff_type t : all) {
      std::ostringstream oss;
      oss << t;
      // operator<< must emit a single whitespace-free token (ToString truncates
      // at the first whitespace, so a space would break the settings round-trip)
      INFO("token=\"" << oss.str() << "\"");
      CHECK(oss.str().find(' ') == std::string::npos);
      std::istringstream iss(oss.str());
      EPA_ff_type back = EPA_ff_type::Test;
      iss >> back;
      CHECK(back == t);
    }
  }

  SECTION("operator>> accepts numeric codes (the SetDefault(size_t(...)) path)")
  {
    EPA_ff_type t = EPA_ff_type::Undefined;
    std::istringstream("0") >> t;
    CHECK(t == EPA_ff_type::point);
    std::istringstream("11") >> t;
    CHECK(t == EPA_ff_type::hcs);
    std::istringstream("14") >> t;
    CHECK(t == EPA_ff_type::WoodSaxon);
  }

  SECTION("operator<< on an out-of-range value is bounded (no infinite recursion)")
  {
    // Pre-fix the fallback printed the value via "<< type", recursing forever.
    std::ostringstream oss;
    oss << static_cast<EPA_ff_type>(12345);
    INFO("out=\"" << oss.str() << "\"");
    CHECK(oss.str().find("12345") != std::string::npos);
    CHECK(oss.str().find("Unknown") != std::string::npos);
  }
}
