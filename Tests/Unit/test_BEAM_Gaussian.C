#include <catch2/catch_all.hpp>

#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "BEAM/Main/Collider_Weight.H"
#include "BEAM/Spectra/Gaussian.H"

#include <cmath>
#include <functional>

using namespace ATOOLS;
using namespace BEAM;

namespace {

  void Boot()
  {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Settings::InitializeMainSettings("");
    if (s_kftable.find(kf_e) == s_kftable.end())
      AddParticle(kf_e, 0.000511, .0, .0, -3, 0, 1, 0, 1, 1, 0, "e-", "e+",
                  "e^{-}", "e^{+}");
    booted = true;
  }

  double Integrate1D(const std::function<double(double)>& f, double a, double b)
  {
    Lambda_Functor functor(&f);
    Gauss_Integrator integrator(&functor);
    return integrator.Integrate(a, b, 1.e-12);
  }

}  // namespace

TEST_CASE("Gaussian beam spectrum is a normalised, truncated Gaussian in x",
          "[BEAM::Gaussian]")
{
  Boot();
  const double energy = 45.5938, spread = 1.e-3, nsigma = 4.;
  Gaussian beam(Flavour(kf_e), energy, spread, nsigma, 0., 1);

  const double x0 = beam.Peak(), xmin = beam.Xmin();

  SECTION("support, peak and nominal energy")
  {
    CHECK(beam.Xmax() == 1.);
    CHECK(x0 == Catch::Approx(1. / (1. + nsigma * spread)));
    CHECK(xmin == Catch::Approx(x0 * (1. - nsigma * spread)));
    CHECK(beam.NominalEnergy() == Catch::Approx(energy).epsilon(1.e-12));
  }

  SECTION("the weight vanishes outside [xmin, 1]")
  {
    CHECK_FALSE(beam.CalculateWeight(xmin * (1. - 1.e-6), 1.));
    CHECK(beam.Weight() == 0.);
    CHECK_FALSE(beam.CalculateWeight(1. + 1.e-6, 1.));
    CHECK(beam.Weight() == 0.);
    CHECK(beam.CalculateWeight(x0, 1.));
    CHECK(beam.Weight() > 0.);
  }

  SECTION("the weight integrates to one over [xmin, 1]")
  {
    auto w = [&](double x) {
      beam.CalculateWeight(x, 1.);
      return beam.Weight();
    };
    CHECK(Integrate1D(w, xmin, 1.) == Catch::Approx(1.).epsilon(1.e-9));
  }

  SECTION("the relative width is the requested spread")
  {
    auto w = [&](double x) {
      beam.CalculateWeight(x, 1.);
      return beam.Weight();
    };
    auto var = [&](double x) { return w(x) * (x - x0) * (x - x0); };
    const double sigma_x = std::sqrt(Integrate1D(var, xmin, 1.));
    // truncation at 4 sigma shrinks the variance by about 0.3%
    CHECK(sigma_x / x0 == Catch::Approx(spread).epsilon(5.e-3));
  }
}

TEST_CASE("Gaussian beam spectrum rejects unphysical parameters",
          "[BEAM::Gaussian]")
{
  Boot();
  const Flavour electron(kf_e);
  CHECK_THROWS(Gaussian(electron, 45., 0., 4., 0., 1));
  CHECK_THROWS(Gaussian(electron, 45., -1.e-3, 4., 0., 1));
  CHECK_THROWS(Gaussian(electron, 45., 1.e-3, 0., 0., 1));
  // nsigma * spread >= 1 would put the lower edge at or below x = 0
  CHECK_THROWS(Gaussian(electron, 45., 0.3, 4., 0., 1));
}

TEST_CASE("Collider_Weight::BoxProbability is the bivariate normal box mass",
          "[BEAM::Collider_Weight]")
{
  SECTION("rho = 0 factorises into two error functions")
  {
    for (double n : {1., 2., 4., 10.}) {
      const double expected = std::erf(n / M_SQRT2) * std::erf(n / M_SQRT2);
      CHECK(Collider_Weight::BoxProbability(0., n, n) ==
            Catch::Approx(expected).epsilon(1.e-9));
    }
    CHECK(Collider_Weight::BoxProbability(0., 2., 3.) ==
          Catch::Approx(std::erf(2. / M_SQRT2) * std::erf(3. / M_SQRT2))
              .epsilon(1.e-9));
  }

  SECTION("reference value against an independent integration")
  {
    // renormalisation erf(n/sqrt2)^2 / P at rho = 0.5, n = 2 is
    // 0.9934118108053424 (scipy, bivariate normal over the box)
    const double n = 2., rho = 0.5;
    const double norm = std::pow(std::erf(n / M_SQRT2), 2) /
                        Collider_Weight::BoxProbability(rho, n, n);
    CHECK(norm == Catch::Approx(0.9934118108053424).epsilon(1.e-8));
  }

  SECTION("symmetric under rho -> -rho, and tends to 1 for a wide box")
  {
    CHECK(Collider_Weight::BoxProbability(0.7, 3., 3.) ==
          Catch::Approx(Collider_Weight::BoxProbability(-0.7, 3., 3.))
              .epsilon(1.e-12));
    CHECK(Collider_Weight::BoxProbability(0.9, 10., 10.) ==
          Catch::Approx(1.).epsilon(1.e-9));
  }
}

TEST_CASE("Correlation weight is normalised on the truncated box",
          "[BEAM::Collider_Weight]")
{
  // The two beams are sampled from their own truncated marginals, so
  // <w> = int int phi(d1) phi(d2) / erf^2 * w(d1, d2) over the box must be 1
  // with w = CorrelationWeight(..., rhonorm).
  for (double rho : {0.5, -0.5, 0.9}) {
    for (double n : {2., 4.}) {
      const double erf2 = std::pow(std::erf(n / M_SQRT2), 2);
      const double rhonorm = erf2 / Collider_Weight::BoxProbability(rho, n, n);
      auto phi = [](double d) {
        return std::exp(-0.5 * d * d) / std::sqrt(2. * M_PI);
      };
      std::function<double(double)> inner = [&](double d1) {
        std::function<double(double)> g = [&](double d2) {
          return phi(d1) * phi(d2) / erf2 *
                 Collider_Weight::CorrelationWeight(d1, d2, rho, rhonorm);
        };
        return Integrate1D(g, -n, n);
      };
      INFO("rho = " << rho << ", n = " << n);
      CHECK(Integrate1D(inner, -n, n) == Catch::Approx(1.).epsilon(1.e-8));
    }
  }

  SECTION("rho = 0 leaves the weight untouched")
  {
    CHECK(Collider_Weight::CorrelationWeight(0.3, -1.2, 0., 1.) ==
          Catch::Approx(1.).epsilon(1.e-14));
  }
}
