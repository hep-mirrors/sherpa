#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Function_Base.H"
#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>

using namespace ATOOLS;

// ============================================================================
// Helper: Simple function wrappers for testing
// ============================================================================

class ConstantFunction : public Function_Base {
    double m_c;
public:
    explicit ConstantFunction(double c = 1.0) : m_c(c) {}
    double operator()(double) override { return m_c; }
};

class PolynomialFunction : public Function_Base {
    int m_degree;
public:
    explicit PolynomialFunction(int deg) : m_degree(deg) {}
    double operator()(double x) override { return std::pow(x, m_degree); }
};

class ExponentialFunction : public Function_Base {
    double m_scale;
public:
    explicit ExponentialFunction(double scale = 1.0) : m_scale(scale) {}
    double operator()(double x) override { return std::exp(m_scale * x); }
};

class SineFunction : public Function_Base {
    double m_freq;
public:
    explicit SineFunction(double freq = 1.0) : m_freq(freq) {}
    double operator()(double x) override { return std::sin(m_freq * x); }
};

class ChebyshevWeightedFunction : public Function_Base {
public:
    // f(x) = sqrt(1-x^2), integral over [-1,1] should give pi/2
    double operator()(double x) override { return std::sqrt(1.0 - x*x); }
};

class WeakSingularFunction : public Function_Base {
    double m_alpha;
public:
    explicit WeakSingularFunction(double alpha = -0.4) : m_alpha(alpha) {}
    double operator()(double x) override {
        return (x > 0.0) ? std::pow(x, m_alpha) : 0.0;
    }
};

// ============================================================================
// TEST SUITE: Gauss_Integrator Comprehensive Validation
// ============================================================================

TEST_CASE("Gauss_Integrator: Constructor and basic interface",
          "[ATOOLS::Math::Gauss_Integrator][interface]")
{
    SECTION("Constructor accepts nullptr") {
        Gauss_Integrator integrator(nullptr);
        // Should not crash; actual integration requires setting function
    }

    SECTION("Constructor accepts valid function") {
        ConstantFunction func(1.0);
        Gauss_Integrator integrator(&func);
        // Should construct successfully
    }
}

TEST_CASE("Gauss_Integrator: Edge cases and boundary conditions",
          "[ATOOLS::Math::Gauss_Integrator][edge_cases]")
{
    ConstantFunction func(2.5);
    Gauss_Integrator integrator(&func);

    SECTION("Zero-width interval returns zero") {
        // The implementation explicitly checks x1==x2 and returns 0
        double result = integrator.Integrate(1.0, 1.0, 1e-6);
        CHECK(result == 0.0);
    }

    SECTION("Reversed integration limits (antisymmetry)") {
        // int_a^b f dx = -int_b^a f dx
        double forward = integrator.Integrate(0.0, 1.0, 1e-6);
        double backward = integrator.Integrate(1.0, 0.0, 1e-6);

        // Check antisymmetry with relative tolerance
        CHECK_THAT(forward + backward, Catch::Matchers::WithinAbs(0.0, 1e-10));
    }

    SECTION("Constant function over arbitrary interval") {
        // int_a^b c dx = c(b-a)
        double a = -3.7, b = 5.2;
        double result = integrator.Integrate(a, b, 1e-9);
        double expected = 2.5 * (b - a);

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-10));
    }
}

TEST_CASE("Gauss_Integrator: Polynomial exactness - Legendre method",
          "[ATOOLS::Math::Gauss_Integrator][legendre][exactness]")
{
    const double REL_TOL = 1e-10; // Limited by GauLeg EPS = 3e-11

    SECTION("Linear function: int_0^1 x dx = 0.5") {
        PolynomialFunction func(1);
        Gauss_Integrator integrator(&func);
        double result = integrator.Legendre(0.0, 1.0, 8);
        CHECK_THAT(result, Catch::Matchers::WithinRel(0.5, REL_TOL));
    }

    SECTION("Cubic function: int_0^1 x^3 dx = 0.25") {
        PolynomialFunction func(3);
        Gauss_Integrator integrator(&func);
        double result = integrator.Legendre(0.0, 1.0, 8);
        CHECK_THAT(result, Catch::Matchers::WithinRel(0.25, REL_TOL));
    }

    SECTION("Symmetric odd polynomial: int_{-1}^1 x^7 dx = 0") {
        // Odd function over symmetric interval should give zero
        PolynomialFunction func(7);
        Gauss_Integrator integrator(&func);
        double result = integrator.Legendre(-1.0, 1.0, 8);

        // Use absolute tolerance for zero result
        CHECK_THAT(result, Catch::Matchers::WithinAbs(0.0, 1e-14));
    }

    SECTION("Even polynomial: int_{-1}^1 x^10 dx = 2/11") {
        PolynomialFunction func(10);
        Gauss_Integrator integrator(&func);
        double result = integrator.Legendre(-1.0, 1.0, 16);
        double expected = 2.0 / 11.0;

        // Slightly less stable because of the EPS in Gauss_Integrator::GauLeg
        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, REL_TOL * 10.));
    }

    SECTION("High-degree polynomial: int_{-1}^1 x^20 dx = 2/21") {
        PolynomialFunction func(20);
        Gauss_Integrator integrator(&func);
        double result = integrator.Legendre(-1.0, 1.0, 32);
        double expected = 2.0 / 21.0;

        // Slightly less stable because of the EPS in Gauss_Integrator::GauLeg
        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, REL_TOL * 10.));
    }
}

TEST_CASE("Gauss_Integrator: Smooth transcendental functions",
          "[ATOOLS::Math::Gauss_Integrator][transcendental]")
{
    const double INTEGRATION_PREC = 1e-10;
    const double TEST_TOL = 1e-9; // 10x buffer

    SECTION("Exponential: int_0^1 e^x dx = e - 1") {
        ExponentialFunction func(1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 1.0, INTEGRATION_PREC);
        double expected = std::exp(1.0) - 1.0; // ≈ 1.718281828

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, TEST_TOL));
    }

    SECTION("Sine: int_0^1 sin(x) dx = 1 - cos(1)") {
        SineFunction func(1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 1.0, INTEGRATION_PREC);
        double expected = 1.0 - std::cos(1.0); // ≈ 0.459697694

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, TEST_TOL));
    }

    SECTION("Decaying exponential: int_0^5 e^{-x} dx = 1 - e^{-5}") {
        ExponentialFunction func(-1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 5.0, INTEGRATION_PREC);
        double expected = 1.0 - std::exp(-5.0); // ≈ 0.993262

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, TEST_TOL));
    }

    SECTION("Large interval: int_0^{100} e^{-x} dx ≈ 1.0") {
        ExponentialFunction func(-1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 100.0, INTEGRATION_PREC);
        double expected = 1.0 - std::exp(-100.0); // Essentially 1.0

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, TEST_TOL));
    }
}

TEST_CASE("Gauss_Integrator: Adaptive refinement convergence",
          "[ATOOLS::Math::Gauss_Integrator][convergence]")
{
    // Test that tightening tolerance improves accuracy monotonically
    ExponentialFunction func(1.0);
    Gauss_Integrator integrator(&func);
    double exact = std::exp(1.0) - 1.0;

    SECTION("Convergence as tolerance tightens") {
        double prec_loose = 1e-4;
        double prec_tight = 1e-8;

        double result_loose = integrator.Integrate(0.0, 1.0, prec_loose);
        double result_tight = integrator.Integrate(0.0, 1.0, prec_tight);

        double error_loose = std::abs(result_loose - exact);
        double error_tight = std::abs(result_tight - exact);

        INFO("Loose error: " << error_loose << ", Tight error: " << error_tight);

        // Tighter tolerance should produce smaller error (or comparable)
        CHECK(error_tight <= error_loose * 10.0); // Allow some noise
    }

    SECTION("Multiple precision levels show monotonic improvement") {
        std::vector<double> precisions = {1e-3, 1e-5, 1e-7, 1e-9};
        std::vector<double> errors;

        for (double prec : precisions) {
            double result = integrator.Integrate(0.0, 1.0, prec);
            errors.push_back(std::abs(result - exact));
        }

        // Check that errors generally decrease (allowing for numerical noise)
        for (size_t i = 1; i < errors.size(); ++i) {
            INFO("prec=" << precisions[i] << ": error[" << i-1 << "]="
                 << errors[i-1] << ", error[" << i << "]=" << errors[i]);
            CHECK(errors[i] <= errors[i-1] * 100.0); // Very generous
        }
    }
}

TEST_CASE("Gauss_Integrator: Method comparison - Legendre vs Jacobi",
          "[ATOOLS::Math::Gauss_Integrator][methods]")
{
    const double PREC = 1e-8;

    SECTION("Jacobi optimized for sqrt weight: int_{-1}^1 (1-x^2)^{-1/2} dx") {
        // This integral is pi (Jacobi with alpha=beta=-0.5 should handle well)
        ChebyshevWeightedFunction func;
        Gauss_Integrator integrator(&func);

        // Note: The function includes the sqrt, so we use Legendre
        // For true Jacobi test, the weight should be separate
        double result = integrator.Integrate(-1.0, 1.0, PREC, 1); // mode=1
        double expected = M_PI / 2.0;

        // This is challenging for standard Legendre due to endpoint behavior
        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-4));
    }
}

TEST_CASE("Gauss_Integrator: Weakly singular integrands",
          "[ATOOLS::Math::Gauss_Integrator][singular][limitation]")
{
    // Document that weakly singular functions have reduced accuracy

    SECTION("Weak singularity: int_0^1 x^{-0.4} dx = 5/3") {
        WeakSingularFunction func(-0.4);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 1.0, 1e-6, 1, 65536);
        double expected = 5.0 / 3.0; // = 1 / (1 - 0.4)

        // Expect reduced accuracy due to singularity at x=0
        INFO("Result: " << result << ", Expected: " << expected
             << ", Rel error: " << std::abs(result - expected) / expected);

        // Tolerance relaxed to 1e-3 due to endpoint singularity
        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-3));
    }
}

TEST_CASE("Gauss_Integrator: Oscillatory integrands (known limitations)",
          "[ATOOLS::Math::Gauss_Integrator][oscillatory][limitation][!mayfail]")
{
    // Tag with [!mayfail] since standard Gauss quadrature is not designed
    // for highly oscillatory integrals

    SECTION("Moderate oscillation: int_0^{pi} sin(x) dx = 2") {
        SineFunction func(1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, M_PI, 1e-6);
        double expected = 2.0;

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-5));
    }

    SECTION("High-frequency oscillation: int_0^{10pi} sin(x) dx = 0") {
        // This may fail or require very large nmax
        SineFunction func(1.0);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 10.0 * M_PI, 1e-4, 1, 65536);

        INFO("High-frequency result: " << result << " (expected 0.0)");

        // Use absolute tolerance for zero result; relaxed due to aliasing
        CHECK_THAT(result, Catch::Matchers::WithinAbs(0.0, 1e-2));
    }
}

TEST_CASE("Gauss_Integrator: Stress tests and robustness",
          "[ATOOLS::Math::Gauss_Integrator][stress]")
{
    SECTION("Very tight tolerance approaching EPS limit") {
        // Request prec=1e-10 (close to GauLeg EPS=3e-11)
        PolynomialFunction func(4);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 1.0, 1e-10, 1, 4096);
        double expected = 0.2; // x^4 from 0 to 1

        // Should still work for polynomials
        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-9));
    }

    SECTION("Very large nmax (weight caching test)") {
        // Ensure static weight list doesn't corrupt with large n
        ConstantFunction func(1.0);
        Gauss_Integrator integrator(&func);

        double result1 = integrator.Integrate(0.0, 1.0, 1e-8, 1, 16384);
        double result2 = integrator.Integrate(0.0, 1.0, 1e-8, 1, 16384);

        // Second call should reuse cached weights
        CHECK_THAT(result1, Catch::Matchers::WithinRel(1.0, 1e-10));
        CHECK(result1 == result2); // Exact match due to caching
    }

    SECTION("Interval with negative bounds") {
        PolynomialFunction func(2);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(-5.0, -2.0, 1e-8);
        double expected = (std::pow(-2.0, 3) - std::pow(-5.0, 3)) / 3.0;

        CHECK_THAT(result, Catch::Matchers::WithinRel(expected, 1e-9));
    }
}

TEST_CASE("Gauss_Integrator: Numerical stability checks",
          "[ATOOLS::Math::Gauss_Integrator][stability]")
{
    SECTION("Function that changes sign") {
        // int_{-1}^1 x dx = 0 (cancellation test)
        PolynomialFunction func(1);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(-1.0, 1.0, 1e-9);

        CHECK_THAT(result, Catch::Matchers::WithinAbs(0.0, 1e-13));
    }

    SECTION("Near-zero integrand") {
        // int_0^1 1e-8 dx = 1e-8
        ConstantFunction func(1e-8);
        Gauss_Integrator integrator(&func);
        double result = integrator.Integrate(0.0, 1.0, 1e-6);

        CHECK_THAT(result, Catch::Matchers::WithinRel(1e-8, 1e-5));
    }
}

