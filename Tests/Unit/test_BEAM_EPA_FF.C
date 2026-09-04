#include "BEAM/Spectra/EPA_FF.H"

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
// Closure test for the b-dependent photon flux.
//
// This is the referee's literal request in point 3 -- "integrating the
// b-dependent flux and comparing it with analytic results" -- carried out
// through the *actual* code path rather than a reimplementation: the
// numerically Fourier-transformed N(x,b) table, its interpolation, the SampleB
// importance weight, and the alpha/pi prefactor that EPA::Initialise applies to
// EPA_FF_Base::N().
//
// EPA_DipoleApprox is the right probe because its form factor is F(Q^2) = 1, so
// EPA_FF_Base::FillTables has to reproduce the elementary transform
//
//     int_0^inf dk_T  k_T^2/(k_T^2 + x^2 m^2)  J_1(b k_T)  =  x m K_1(x m b)
//
// exactly. Everything the referee questioned about the flux machinery -- the
// oscillatory integrator, the kernel, the table, the sampling weight and the
// prefactors -- is thereby tested against a closed form, and the test itself
// assumes nothing about the point-like limit.
//
// Reference results (paper, Sec. 2.2):
//   n_PL(x,b) = (Z^2 alpha/pi^2) x m^2
//               [ K_1^2(chi) + K_0^2(chi)/gamma^2 ] .
// In the conventional ultrarelativistic limit, gamma -> infinity, this gives
//   n_PL,int(x) = int_{b_min}^{inf} d^2b n_PL
//               = (2 Z^2 alpha/pi x) [ chi K_0 K_1 - (chi^2/2)(K_1^2 - K_0^2) ] .
// EPA_IonApproxIntegrated implements that transverse-only limit. The finite-
// gamma closure below keeps the K_0 term over precisely the range where the
// implementation does. The indefinite-integral form accounts exactly for the
// finite upper limit imposed by SampleB's b_max/chi_max cap.
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

  // EPA settings the form-factor constructors read. Mirrors the relevant part
  // of EPA::RegisterDefaults, at the production defaults; the two
  // result-directory keys are pinned so that EPA_FF_Base leaves m_respath empty
  // and no table cache is read or written.
  constexpr size_t s_xbins = 200, s_bbins = 100;
  constexpr double s_xmin = 1.e-5, s_xmax = 1.;
  constexpr double s_bmin = 0.3, s_bmax = 1.e3, s_bthr = 10., s_chimax = 100.;
  constexpr double s_radius = 0.88; // proton radius in fm, cf. Tab. 2
  constexpr double s_mass = 0.938272;
  constexpr double s_alpha = 1. / 137.03599976;
  constexpr double s_ebeam = 6500.; // LHC proton beam, sets gamma

  void Boot()
  {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    ATOOLS::Settings::InitializeMainSettings("");
    if (ATOOLS::s_kftable.find(kf_p_plus) == ATOOLS::s_kftable.end())
      ATOOLS::AddParticle(kf_p_plus, s_mass, s_radius, 0., 3, 1, true, 1, "P+",
                          "P+");
    // The kf table is shared across the whole test binary and Catch2 randomises
    // test order, so another test case may already have registered the proton --
    // the REMNANTS tests do, with radius 0 and mass 0.938. AddParticle above is
    // then skipped, so pin both properties rather than relying on who ran first:
    // EPA_FF_Base divides by the radius, and every reference value below is
    // built from s_mass, so a 3e-4 mass mismatch would quietly bias the whole
    // comparison. Pinning is safe in the other direction too: Beam_Base's
    // constructor is the only other reader of the radius in the suite (via the
    // REMNANTS Test_Beam classes) and merely stores it, and the REMNANTS tests
    // assert on neither.
    ATOOLS::Flavour(kf_p_plus).SetMass(s_mass);
    ATOOLS::Flavour(kf_p_plus).SetRadius(s_radius);
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
  // SECTION, and filling the 201x101 table costs ~2*10^4 oscillatory integrals.
  // Leaking also avoids ordering the destructor against rpa/Settings teardown.
  EPA_DipoleApprox& Flux()
  {
    static EPA_DipoleApprox* ff = nullptr;
    if (!ff) {
      Boot();
      // The beam energy is a constructor argument, as it has to be: the
      // Lorentz factor must be known before a derived constructor fills its
      // N(x,b) table, which is what EPA_IonApprox does. The first section below
      // asserts it arrived.
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

  // Transverse term of Eq. (11) in the normalisation FillTables stores, i.e.
  // including the d^2b = 2 pi b db Jacobian and without the alpha/pi prefactor.
  // The table carries no longitudinal term: it descends from C(Q^2) = 0.
  double AnalyticFlux(double x, double b_phys)
  {
    const double chi = x * s_mass * b_phys;
    return 2. * b_phys * x * s_mass * s_mass * ATOOLS::sqr(K1(chi));
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

  SECTION("the beam energy reaches gamma through the constructor")
  {
    // Regression guard for the ordering: gamma has to be settled by the end of
    // EPA_FF_Base's constructor, because EPA_IonApprox fills a gamma-dependent
    // table inside its own. There is deliberately no setter to install it late.
    CHECK(Flux().Gamma() == Catch::Approx(gamma));
  }

  SECTION("the numerical transform reproduces x m K_1(x m b) at table nodes")
  {
    // On a node the only error left is the oscillatory integrator's own
    // tolerance (Bessel_Integrator default 1e-3), so this isolates the
    // transform from the interpolation tested below.
    for (size_t i : {40u, 80u, 120u}) { // x = 1e-4, 1e-3, 1e-2
      const double x = xaxis.x(i);
      for (size_t j : {10u, 30u, 50u, 70u, 90u}) {
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

  SECTION("interpolation between nodes stays at the per-mille level")
  {
    // Probed at the *geometric mid-bin* in both axes, i.e. weight 1/2, where
    // the linear-interpolation error w(1-w) peaks. Sampling at fixed x values
    // instead would make the result depend on how close those happen to fall
    // to a node, which changes with the bin count and would misreport the
    // effect of changing it.
    //
    // The residual here is the convexity of N across a bin: the *metric* half
    // of the error was removed by giving axis::weight the log-axis fraction,
    // and what is left falls as O(h^2) in the node spacing -- which is why the
    // default xBins was raised to 200. Worst case measured at that default:
    // +4.6e-4, in the steepest corner of the grid (largest x, largest b).
    const auto mid = [](const ATOOLS::axis& a, size_t i) {
      return std::sqrt(a.x(i) * a.x(i + 1));
    };
    for (size_t i : {40u, 80u, 120u}) { // mid-bin above x = 1e-4, 1e-3, 1e-2
      const double x = mid(xaxis, i);
      for (size_t j : {10u, 30u, 50u, 70u, 90u}) {
        const double b_dimless = mid(baxis, j) / R;
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

  SECTION("int d^2b of the sampled flux reproduces n_PL,int(x)")
  {
    // Midpoint quadrature over ran: SampleB returns db/dran, so the mean of
    // N(x,ran) over a uniform ran is exactly int db (dn/db) over the sampled
    // range. Multiplying by alpha/pi is what EPA::CalculateWeight does.
    const int N = 20000;
    for (double x : {1.e-4, 3.e-4, 1.e-3, 3.e-3, 1.e-2, 3.e-2}) {
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
      const double want = pref * 2. / x *
                          ((A1(chi_b) - A1(chi_a)) +
                           (A0(chi_b) - A0(chi_t)) / ATOOLS::sqr(gamma));

      INFO("x=" << x << " chi in [" << chi_a << ", " << chi_b
                << "] got=" << got << " want=" << want
                << " rel=" << (got / want - 1.));
      CHECK_THAT(got, Catch::Matchers::WithinRel(want, 3.e-3));
    }
  }
}
