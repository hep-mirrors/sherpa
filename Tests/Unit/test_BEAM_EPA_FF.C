#include "BEAM/Spectra/EPA_FF.H"

#include "Test_Setup.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

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
      EPA_ff_type::protonSachs,  EPA_ff_type::Test,
      EPA_ff_type::Undefined};

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
    std::istringstream("17") >> t;
    CHECK(t == EPA_ff_type::protonSachs);
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

////////////////////////////////////////////////////////////////////////////////
// Closure test for the b-dependent photon flux, through the production code
// path: the numerically Fourier-transformed N(x,b) table, its interpolation,
// the SampleB importance weight and the alpha/pi prefactor that EPA applies to
// EPA_FF_Base::N().
//
// EPA_DipoleApprox has the form factor F(Q^2) = 1, so EPA_FF_Base::FillTables
// has to reproduce the elementary transform
//
//     int_0^inf dk_T  k_T^2/(k_T^2 + x^2 m^2)  J_1(b k_T)  =  x m K_1(x m b) ,
//
// times (1-x) for the exact virtuality. Its b-integral is compared with the
// point-like flux of a charge Z (Bertulani:1987tz),
//   n_PL(x,b) = (Z^2 alpha/pi^2) x m^2 [ K_1^2(chi) + K_0^2(chi)/gamma^2 ] ,
// chi = x m b, which for gamma -> infinity integrates from b_min to infinity to
//   n_PL,int(x) = (2 Z^2 alpha/pi x) [ chi K_0 K_1 - (chi^2/2)(K_1^2 - K_0^2) ]
// as in EPA_IonApproxIntegrated. The closure below keeps the K_0^2/gamma^2
// term where the implementation does, beyond the point-like threshold, and
// integrates up to the upper limit of SampleB.
////////////////////////////////////////////////////////////////////////////////

namespace {

  double K0(double chi) { return ATOOLS::SF.Kn(0, chi); }
  double K1(double chi) { return ATOOLS::SF.Kn(1, chi); }

  // Antiderivatives, verified by differentiation using K0' = -K1 and
  // K1' = -K0 - K1/chi:
  //   d/dchi A1 = chi K_1^2 ,   d/dchi A0 = chi K_0^2 .
  // A1(chi -> inf) = 0, so -A1(chi) reproduces the bracket of n_PL,int above.
  double A1(double chi)
  {
    const double k0 = K0(chi), k1 = K1(chi);
    return 0.5 * chi * chi * (k1 * k1 - k0 * k0) - chi * k0 * k1;
  }
  double A0(double chi)
  {
    const double k0 = K0(chi), k1 = K1(chi);
    return 0.5 * chi * chi * (k0 * k0 - k1 * k1);
  }

  // EPA settings the form-factor constructors read, as EPA::RegisterDefaults
  // registers them except for a coarser table (production: 200 x 100 bins) that
  // keeps the test fast; 4 x-bins per decade put x = 1e-4, 1e-3, 1e-2 and
  // 10^-0.5 on nodes. The two result-directory keys are pinned so that
  // EPA_FF_Base leaves m_respath empty and no table cache is read or written.
  constexpr size_t s_xbins = 20, s_bbins = 50;
  constexpr double s_xmin = 1.e-5, s_xmax = 1.;
  constexpr double s_bmin = 0.3, s_bmax = 1.e3, s_bthr = 10., s_chimax = 100.;
  constexpr double s_mass = Test_Setup::proton_mass;
  constexpr double s_alpha = 1. / 137.03599976;
  constexpr double s_ebeam = 6500.; // LHC proton beam, sets gamma

  void Boot()
  {
    static bool booted = false;
    if (booted) return;
    Test_Setup::InitGlobals();
    ATOOLS::Settings::InitializeMainSettings("");
    Test_Setup::RegisterParticles();
    auto& s = ATOOLS::Settings::GetMainSettings();
    s["GENERATE_RESULT_DIRECTORY"].SetDefault(false);
    s["RESULT_DIRECTORY"].SetDefault("");
    auto epa = s["EPA"];
    epa["Q2Max"].SetDefault(1.);
    epa["Q2Min"].SetDefault(-1.);
    epa["xMin"].SetDefault(s_xmin);
    epa["xMax"].SetDefault(s_xmax);
    epa["xBins"].SetDefault(static_cast<int>(s_xbins));
    epa["bBins"].SetDefault(static_cast<int>(s_bbins));
    epa["bMin"].SetDefault(s_bmin);
    epa["bMax"].SetDefault(s_bmax);
    epa["bThreshold"].SetDefault(s_bthr);
    epa["chiMax"].SetDefault(s_chimax);
    booted = true;
  }

  // Built once and deliberately leaked: Catch2 re-runs a TEST_CASE body per
  // SECTION, and each table node costs an oscillatory integral. Leaking also
  // avoids ordering the destructor against rpa/Settings teardown.
  EPA_DipoleApprox& Flux()
  {
    static EPA_DipoleApprox* ff = nullptr;
    if (!ff) {
      Boot();
      ff = new EPA_DipoleApprox(ATOOLS::Flavour(kf_p_plus), 1, s_ebeam);
    }
    return *ff;
  }

  double ProtonR() // beam radius in 1/GeV, as the form factor actually holds it
  {
    return Flux().Radius();
  }

  // Upper sampling limit in units of R, as in EPA_FF_Base::SampleB.
  double bmaxx_of(double x, double R)
  {
    return std::max(s_bmin, std::min(s_bmax, s_chimax / (x * s_mass * R)));
  }

  // Invert b(ran) of SampleB, so a chosen impact parameter can be probed
  // through the public interface.
  double ran_of_b(double b_dimless, double x, double R)
  {
    const double norm = 0.5 * std::log((ATOOLS::sqr(bmaxx_of(x, R)) + 1.) /
                                       (ATOOLS::sqr(s_bmin) + 1.));
    return 0.5 *
           std::log((ATOOLS::sqr(b_dimless) + 1.) / (ATOOLS::sqr(s_bmin) + 1.)) /
           norm;
  }

  // Read the stored N(x,b) back through the public interface: N() returns
  // table * weight and SampleB returns that same weight for the same ran.
  double TabulatedFlux(double x, double b_dimless, double R, double& b_phys)
  {
    const double ran = ran_of_b(b_dimless, x, R);
    const double wt = EPA_FF_Base::SampleB(x, ran, s_bmin, s_bmax, s_chimax,
                                           s_mass, R, b_phys);
    return Flux().N(x, ran) / wt;
  }

  // Transverse term of n_PL in the normalisation FillTables stores, i.e.
  // including the d^2b = 2 pi b db Jacobian and without the alpha/pi prefactor.
  // The table carries no longitudinal term: it descends from C(Q^2) = 0.
  // With the exact virtuality the transform evaluates to (1-x) x m K_1(chi),
  // and chi = x m b is independent of that choice.
  double AnalyticFlux(double x, double b_phys)
  {
    const double chi = x * s_mass * b_phys;
    return 2. * b_phys * x * s_mass * s_mass * (1. - x) * ATOOLS::sqr(K1(chi));
  }

} // namespace

TEST_CASE("EPA b-dependent flux closes against the analytic point-like result",
          "[BEAM::EPA_FF][closure]")
{
  const double R = ProtonR();
  const double gamma = s_ebeam / s_mass;
  const double pref = s_alpha / M_PI;

  // Node positions of the two table axes, constructed exactly as
  // EPA_FF_Base::FillTables does (the b axis stops at the point-like
  // threshold, min(b_thr, b_max) * R = 10 R here).
  const ATOOLS::axis xaxis(s_xbins, s_xmin, s_xmax, ATOOLS::axis_mode::log);
  const ATOOLS::axis baxis(s_bbins, s_bmin * R, s_bthr * R,
                           ATOOLS::axis_mode::log);

  // the beam energy is a constructor argument, so that a derived constructor
  // can fill a gamma-dependent table, as EPA_IonApprox does
  CHECK(Flux().Gamma() == Catch::Approx(gamma));

  SECTION("the numerical transform reproduces x m K_1(x m b) at table nodes")
  {
    // on a node the only error is the tolerance of the oscillatory integrator
    // (Bessel_Integrator default 1e-3)
    for (size_t i : {4u, 8u, 12u}) { // x = 1e-4, 1e-3, 1e-2
      const double x = xaxis.x(i);
      for (size_t j : {5u, 15u, 25u, 35u, 45u}) {
        const double b_dimless = baxis.x(j) / R;
        double b_phys = -1.;
        const double got = TabulatedFlux(x, b_dimless, R, b_phys);
        const double want = AnalyticFlux(x, b_phys);
        INFO("x=" << x << " b/R=" << b_dimless << " chi=" << x * s_mass * b_phys
                  << " got=" << got << " want=" << want
                  << " rel=" << (got / want - 1.));
        CHECK_THAT(got, Catch::Matchers::WithinRel(want, 2.e-3));
      }
    }
  }

  SECTION("the (1-x) is present, probed where it is large")
  {
    // node 18 is x = 10^-0.5, where the factor is 0.684. Only the inner b nodes,
    // chi <~ 2.3, are probed: further out the oscillatory integration loses the
    // exponentially small transform in the cancellation between its partial
    // sums. Production stays below chi = 1.
    const double x = xaxis.x(18);
    REQUIRE(x == Catch::Approx(0.31622776601683794));
    for (size_t j : {5u, 15u, 25u}) {
      const double b_dimless = baxis.x(j) / R;
      double b_phys = -1.;
      const double got = TabulatedFlux(x, b_dimless, R, b_phys);
      const double chi = x * s_mass * b_phys;
      // Spelled out rather than taken from AnalyticFlux, so that dropping the
      // factor in the code cannot be masked by dropping it in the reference.
      const double want_exact =
          2. * b_phys * x * s_mass * s_mass * (1. - x) * ATOOLS::sqr(K1(chi));
      INFO("x=" << x << " b/R=" << b_dimless << " chi=" << chi << " got=" << got
                << " want=" << want_exact << " rel=" << (got / want_exact - 1.)
                << " (without the (1-x): " << want_exact / (1. - x) << ")");
      CHECK_THAT(got, Catch::Matchers::WithinRel(want_exact, 2.e-3));
    }
  }

  SECTION("int d^2b of the sampled flux reproduces n_PL,int(x)")
  {
    // Midpoint quadrature over ran: SampleB returns db/dran, so the mean of
    // N(x,ran) over a uniform ran is exactly int db (dn/db) over the sampled
    // range. Multiplying by alpha/pi is what EPA::CalculateWeight does. x sits
    // on table nodes, from 1e-4 to 10^-1.5 in half decades.
    const int N = 20000;
    for (size_t i : {4u, 6u, 8u, 10u, 12u, 14u}) {
      const double x = xaxis.x(i);
      double sum = 0.;
      for (int i = 0; i < N; ++i) sum += Flux().N(x, (i + 0.5) / N);
      const double got = pref * sum / N;

      const double bmaxx = bmaxx_of(x, R);
      const double chi_a = x * s_mass * s_bmin * R;
      const double chi_b = x * s_mass * bmaxx * R;
      // Above b_thr = 10 R the flux switches to the closed-form point-like
      // branch, which adds the longitudinal K_0^2/gamma^2 term; below it the
      // table carries the transverse term only.
      const double chi_t = std::min(chi_b, x * s_mass * s_bthr * R);
      // (1-x) as in AnalyticFlux, b-independent and hence outside the integral.
      const double want = pref * 2. * (1. - x) / x *
                          ((A1(chi_b) - A1(chi_a)) +
                           (A0(chi_b) - A0(chi_t)) / ATOOLS::sqr(gamma));

      INFO("x=" << x << " chi in [" << chi_a << ", " << chi_b
                << "] got=" << got << " want=" << want
                << " rel=" << (got / want - 1.));
      CHECK_THAT(got, Catch::Matchers::WithinRel(want, 3.e-3));
    }
  }
}
