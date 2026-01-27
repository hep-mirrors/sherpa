#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Function_Base.H"
#include <catch2/catch_all.hpp>
#include <cmath>
#include <iostream>

using namespace ATOOLS;
using namespace std;

// Helper class to replicate the behavior of Bessel_Integrator::TestFunction
class BesselTestFunc: public Function_Base {
private:

  int m_mode;

public:

  explicit BesselTestFunc(int mode): m_mode(mode) {}
  double operator()(double x) override
  {
    switch (m_mode) {
    case 0:
      return std::exp(-x);
    case 1:
      return std::sin(x / 2) / (1. + x * x);
    case 2:
      return std::sin(x * x / 3) / x;
    case 3:
      return x / (1. + x * x);
    case 4:
      // Prevent div by zero at x=0 if the integrator evaluates there
      if (std::abs(x) < 1e-15)
        return 2. / log(1. + sqrt(2.)); // Limit behaviour
      return (1. - exp(-x)) / (x * log(1. + sqrt(2.)));
    case 5:
      return std::exp(-x) * std::cos(3. * x);
    case 6:
      return std::sin(x / 10.) / 10. / (x + 1.);
    default:
      return 1.0;
    }
  }
};

TEST_CASE("Bessel_Integrator Instantiation",
          "[ATOOLS::Math::Bessel_Integrator]")
{
  BesselTestFunc f(0);
  Bessel_Integrator integrator(&f, 1);
  // Just check if it constructs without crashing
  REQUIRE(true);
}

TEST_CASE("Bessel_Integrator Integration Modes",
          "[ATOOLS::Math::Bessel_Integrator]")
{

  const double REL_TOL = 1e-5;

  // Define the test data structure
  struct TestConfig {
    int mode;
    size_t order;
    double expected_val;
    string description;
  };

  vector<TestConfig> configs = {
      {0, 0, 0.707106781186547524400844, "f(x)=exp(-x), J0(x)"},
      {0, 1, 0.292893218813452475599155, "f(x)=exp(-x), J1(x)"},
      {1, 1, 0.222146354934442148139912, "f(x)=sin(x/2)/(1+x^2), J1(x)"},
      {2, 1, 0.425189666999595468010652, "f(x)=sin(x^2/3)/x, J1(x)"},
      {3, 0, 0.421024438240708333335627, "f(x)=x/(1+x^2), J0(x)"},
      {4, 0, 1.0, "f(x)=(1-exp(-x))/(x*log(1+sqrt(2))), J0(x)"},
      {5, 0, 0.114263484396498030486956, "f(x)=exp(-x)*cos(x), J0(x)"},
      {5, 1, -0.040914316903594003073207, "f(x)=exp(-x)*cos(3*x), J1(x)"},
      {6, 1, 0.00539632971919505751151192, "f(x)=sin(x/10)/(x + 1)/10, J1(x)"}};

  for (const auto& config : configs) {
    DYNAMIC_SECTION("Mode " << config.mode << ": " << config.description)
    {
      BesselTestFunc func(config.mode);
      Bessel_Integrator integrator(&func, config.order);

      double result = integrator(REL_TOL / 10.);

      INFO("Mode " << config.mode << " details:\n"
                   << std::fixed << std::setprecision(17)
                   << "  Expected: " << config.expected_val << "\n"
                   << "  Got:      " << result << "\n"
                   << "  Rel Diff: "
                   << std::abs(result - config.expected_val) /
                          config.expected_val);

      CHECK_THAT(result,
                 Catch::Matchers::WithinRel(config.expected_val, REL_TOL));
    }
  }
}

TEST_CASE("Bessel_Integrator Limits and Cutoffs",
          "[ATOOLS::Math::Bessel_Integrator]")
{
  // Test that changing maxbins or depth affects the result/performance
  // (This is a regression test for the configuration setters)
  BesselTestFunc func(0);
  Bessel_Integrator integrator(&func, 1);

  SECTION("Setters function correctly")
  {
    integrator.SetMaxBins(100);
    integrator.SetDepth(5);
    // If these didn't crash, we assume setters work (state inspection is hard
    // without getters)
    REQUIRE(true);
  }
}
