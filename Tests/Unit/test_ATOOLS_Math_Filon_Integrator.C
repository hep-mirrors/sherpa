#include "ATOOLS/Math/Filon_Integrator.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/MathTools.H"

#include <catch2/catch_all.hpp>
#include <cmath>
#include <functional>

using namespace ATOOLS;

// The Filon quadrature must be accurate even when sin(k r) oscillates many
// times across the grid, as long as the grid resolves the smooth factor.
TEST_CASE("Filon_Integrator matches closed-form transforms", "[ATOOLS::Filon]")
{
  SECTION("exp(-r): int_0^inf e^{-r} sin(k r) dr = k/(1+k^2)")
  {
    std::function<double(double)> f = [](double r) { return std::exp(-r); };
    Lambda_Functor fl(&f);
    Filon_Integrator filon(&fl, 0., 40., 2000); // tail e^{-40} negligible
    for (double k : {0.25, 1.0, 5.0, 20.0, 100.0}) {
      const double exact = k / (1. + k * k);
      INFO("k = " << k << " exact = " << exact);
      CHECK_THAT(filon.SineTransform(k),
                 Catch::Matchers::WithinAbs(exact, 1e-9));
    }
  }

  SECTION("exp(-r): int_0^inf e^{-r} cos(k r) dr = 1/(1+k^2)")
  {
    std::function<double(double)> f = [](double r) { return std::exp(-r); };
    Lambda_Functor fl(&f);
    Filon_Integrator filon(&fl, 0., 40., 2000);
    for (double k : {0.25, 1.0, 5.0, 20.0, 100.0}) {
      const double exact = 1. / (1. + k * k);
      INFO("k = " << k << " exact = " << exact);
      CHECK_THAT(filon.CosineTransform(k),
                 Catch::Matchers::WithinAbs(exact, 1e-9));
    }
  }

  SECTION("r*exp(-r): int_0^inf r e^{-r} sin(k r) dr = 2*k/(1+k^2)^2")
  {
    std::function<double(double)> f = [](double r) { return r*std::exp(-r); };
    Lambda_Functor fl(&f);
    Filon_Integrator filon(&fl, 0., 40., 2000);
    for (double k : {0.25, 1.0, 5.0, 20.0, 100.0}) {
      const double exact = 2. * k / std::pow(1. + k * k, 2);
      INFO("k = " << k << " exact = " << exact);
      CHECK_THAT(filon.SineTransform(k),
                 Catch::Matchers::WithinAbs(exact, 1e-9));
    }
  }

  SECTION("constant f: int_0^L sin(k r) dr = (1 - cos(k L))/k")
  {
    std::function<double(double)> f = [](double) { return 1.; };
    Lambda_Functor fl(&f);
    const double L = 3.;
    Filon_Integrator filon(&fl, 0., L, 100);
    for (double k : {0.5, 2.0, 7.0, 30.0}) {
      const double exact = (1. - std::cos(k * L)) / k;
      INFO("k = " << k << " exact = " << exact);
      CHECK_THAT(filon.SineTransform(k),
                 Catch::Matchers::WithinAbs(exact, 1e-12));
    }
  }

  SECTION("CosineTransform(0) equals the plain integral")
  {
    std::function<double(double)> f = [](double r) { return std::exp(-r); };
    Lambda_Functor fl(&f);
    Filon_Integrator filon(&fl, 0., 40., 2000);
    CHECK_THAT(filon.CosineTransform(0.),
               Catch::Matchers::WithinAbs(1. - std::exp(-40.), 1e-9));
  }

  SECTION("Woods-Saxon-like integrand: coarse grid already accurate")
  {
    // g(r) = r * rho(r) with a Pb-like Fermi profile (units of 1/GeV).
    const double R = 33., d = 2.8;
    std::function<double(double)> g = [&](double r) {
      return r / (1. + std::exp((r - R) / d));
    };
    Lambda_Functor gl(&g);
    const double rmax = R + 16. * d;
    Filon_Integrator coarse(&gl, 0., rmax, 1024);
    Filon_Integrator fine(&gl, 0., rmax, 8192);
    // Both grids resolve rho(r); agreement => the default 1024-interval grid is
    // already converged, including where sin(q r) oscillates fast.
    for (double q : {0.5, 1.0, 5.0, 10.0}) {
      INFO("q = " << q << " coarse = " << coarse.SineTransform(q)
                  << " fine = " << fine.SineTransform(q));
      CHECK_THAT(coarse.SineTransform(q),
                 Catch::Matchers::WithinRel(fine.SineTransform(q), 1e-6));
    }
  }
}
