#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "ATOOLS/Math/Special_Functions.H"

using namespace ATOOLS;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("Bessel J0: Small arguments (series regime)", "[bessel][j0]") {
    Special_Functions SF;

    SECTION("x = 0 (exact value)") {
        REQUIRE_THAT(SF.J0(0.0), WithinAbs(1.0, 1e-15));
    }

    SECTION("x \to 0 (series convergence)") {
        // Reference values from Wolfram Alpha / NIST DLMF
        REQUIRE_THAT(SF.J0(1e-12), WithinRel(1.0, 1e-14));
        REQUIRE_THAT(SF.J0(1e-6), WithinRel(0.9999999999997500, 1e-12));
    }

    SECTION("Small x < 1 (series regime)") {
        // J_0(0.1) = 0.997501562066040...
        REQUIRE_THAT(SF.J0(0.1), WithinRel(0.997501562066040, 1e-8));

        // J_0(0.5) = 0.938469807240813...
        REQUIRE_THAT(SF.J0(0.5), WithinRel(0.938469807240813, 1e-9));
    }

    SECTION("Intermediate 1 < x < 8 (rational approximation)") {
        // J_0(1.0) = 0.765197686557967...
        REQUIRE_THAT(SF.J0(1.0), WithinRel(0.765197686557967, 1e-8));

        // J_0(2.0) = 0.223890779141236...
        REQUIRE_THAT(SF.J0(2.0), WithinRel(0.223890779141236, 1e-7));

        // Near first zero: J_0(2.404825557695773) \approx 0
        // Absolute tolerance required near zeros
        REQUIRE_THAT(SF.J0(2.404825557695773), WithinAbs(0.0, 1e-8));

        // J_0(5.0) = -0.17759677131434...
        REQUIRE_THAT(SF.J0(5.0), WithinRel(-0.17759677131434, 1e-7));

        // J_0(7.5) = -0.2663396578804...
        REQUIRE_THAT(SF.J0(7.5), WithinRel(0.2663396578804, 1e-8));
    }
}

TEST_CASE("Bessel J0: Large arguments (asymptotic regime)", "[bessel][j0]") {
    Special_Functions SF;

    SECTION("Transition region x \\approx 8") {
        // J_0(8.0) = 0.1716508071375539...
        REQUIRE_THAT(SF.J0(8.0), WithinRel(0.1716508071375539, 1e-9));

        // J_0(10.0) = -0.2459357644513565...
        REQUIRE_THAT(SF.J0(10.0), WithinRel(-0.2459357644513565, 1e-9));
    }

    SECTION("Large x (asymptotic expansion)") {
        // J_0(20.0) = 0.167024664340583...
        REQUIRE_THAT(SF.J0(20.0), WithinRel(0.167024664340583, 1e-9)); // fails

        // J_0(50.0) = -0.055812327669252...
        REQUIRE_THAT(SF.J0(50.0), WithinRel(0.055812327669252, 1e-10));

        // J_0(100.0) = 0.019985850304223...
        REQUIRE_THAT(SF.J0(100.0), WithinRel(0.019985850304223, 5e-10));
    }

    SECTION("Very large x (asymptotic accuracy degrades)") {
        // J_0(200.0) = -0.01543743993057...
        // Looser tolerance due to sin/cos error propagation
        REQUIRE_THAT(SF.J0(200.0), WithinRel(-0.01543743993057, 1e-8));
    }
}

TEST_CASE("Bessel J0: Near zeros", "[bessel][j0][zeros]") {
    Special_Functions SF;

    // First few zeros of J_0 (from NIST DLMF)
    // Use absolute tolerance near zeros - relative error is unbounded

    SECTION("First zero") {
        double z1 = 2.404825557695773;
        REQUIRE_THAT(SF.J0(z1), WithinAbs(0.0, 1e-8));
    }

    SECTION("Second zero") {
        double z2 = 5.520078110286311;
        REQUIRE_THAT(SF.J0(z2), WithinAbs(0.0, 1e-8));
    }

    SECTION("Third zero") {
        double z3 = 8.653727912911012;
        REQUIRE_THAT(SF.J0(z3), WithinAbs(0.0, 1e-8));
    }

    SECTION("Large zero") {
        double z10 = 30.63460646843198;
        REQUIRE_THAT(SF.J0(z10), WithinAbs(0.0, 1e-7));
    }
}

TEST_CASE("Bessel J0: Negative arguments", "[bessel][j0][negative]") {
    Special_Functions SF;

    // J_0 is an even function: J_0(-x) = J_0(x)
    SECTION("Even function property") {
        REQUIRE_THAT(SF.J0(-1.0), WithinRel(SF.J0(1.0), 1e-15));
        REQUIRE_THAT(SF.J0(-5.5), WithinRel(SF.J0(5.5), 1e-14));
        REQUIRE_THAT(SF.J0(-20.0), WithinRel(SF.J0(20.0), 1e-12));
    }
}

TEST_CASE("Bessel J1: Small arguments", "[bessel][j1]") {
    Special_Functions SF;

    SECTION("x = 0 (exact value)") {
        REQUIRE_THAT(SF.J1(0.0), WithinAbs(0.0, 1e-16));
    }

    SECTION("x \\to 0 behavior: J_1(x) ~ x/2") {
        REQUIRE_THAT(SF.J1(1e-12), WithinAbs(5e-13, 1e-12));
        REQUIRE_THAT(SF.J1(1e-6), WithinRel(4.999999999999375e-7, 1e-12)); // fails
    }

    SECTION("Small x (series regime)") {
        // J_1(0.1) = 0.04993752603624200...
        REQUIRE_THAT(SF.J1(0.1), WithinRel(0.04993752603624200, 1e-9)); // fails

        // J_1(0.5) = 0.2422684576748739...
        REQUIRE_THAT(SF.J1(0.5), WithinRel(0.2422684576748739, 1e-10));

        // J_1(1.0) = 0.4400505857449335...
        REQUIRE_THAT(SF.J1(1.0), WithinRel(0.4400505857449335, 1e-9));
    }
}

TEST_CASE("Bessel J1: Intermediate arguments", "[bessel][j1]") {
    Special_Functions SF;

    SECTION("Transition region") {
        // J_1(2.0) = 0.5767248077568734...
        REQUIRE_THAT(SF.J1(2.0), WithinRel(0.5767248077568734, 1e-9)); // fails

        // J_1(5.0) = -0.3275791375914652...
        REQUIRE_THAT(SF.J1(5.0), WithinRel(-0.3275791375914652, 1e-8));

        // J_1(10.0) = 0.04347274616886144...
        REQUIRE_THAT(SF.J1(10.0), WithinRel(0.04347274616886144, 1e-9));
    }
}

TEST_CASE("Bessel J1: Large arguments", "[bessel][j1]") {
    Special_Functions SF;

    SECTION("Large x ≥ 100 (asymptotic regime)") {
        // J_1(100.0) = -0.07714535201411216...
        REQUIRE_THAT(SF.J1(100.0), WithinRel(-0.07714535201411216, 5e-9));

        // J_1(200.0) = 0.05430453818237822...
        REQUIRE_THAT(SF.J1(200.0), WithinRel(-0.05430453818237822, 1e-8)); // fails
    }

    SECTION("Very large x (reduced accuracy expected)") {
        // J_1(500.0) = 0.01047261347037229...
        // Looser tolerance due to implementation limitations
        REQUIRE_THAT(SF.J1(500.0), WithinRel(0.01047261347037229, 1e-7));
    }
}

TEST_CASE("Bessel J1: Near zeros", "[bessel][j1][zeros]") {
    Special_Functions SF;

    SECTION("First zero") {
        double z1 = 3.831705970207515;
        REQUIRE_THAT(SF.J1(z1), WithinAbs(0.0, 1e-8));
    }

    SECTION("Second zero") {
        double z2 = 7.015586669815619;
        REQUIRE_THAT(SF.J1(z2), WithinAbs(0.0, 1e-8));
    }

    SECTION("Third zero") {
        double z3 = 10.17346813506272;
        REQUIRE_THAT(SF.J1(z3), WithinAbs(0.0, 1e-8));
    }
}

TEST_CASE("Bessel J1: Negative arguments", "[bessel][j1][negative]") {
    Special_Functions SF;

    // J_1 is an odd function: J_1(-x) = -J_1(x)
    SECTION("Odd function property") {
        REQUIRE_THAT(SF.J1(-1.0), WithinRel(-SF.J1(1.0), 1e-15));
        REQUIRE_THAT(SF.J1(-5.5), WithinRel(-SF.J1(5.5), 1e-14));
        REQUIRE_THAT(SF.J1(-100.0), WithinRel(-SF.J1(100.0), 1e-12));
    }
}

TEST_CASE("Bessel functions: Stress tests", "[bessel][stress]") {
    Special_Functions SF;

    SECTION("J_0: Oscillatory behavior preservation") {
        // Check multiple periods to verify asymptotic formula consistency
        std::vector<double> test_points = {20.0, 40.0, 60.0, 80.0, 100.0};

        for (double x : test_points) {
            double result = SF.J0(x);
            // Should be bounded: |J_0(x)| < 1
            REQUIRE(std::abs(result) < 1.0);
            // Should decay roughly as 1/sqrt(x)
            REQUIRE(std::abs(result) < 1.0 / std::sqrt(x) + 0.1);
        }
    }

    SECTION("J_1: Monotonicity checks near extrema") {
        // J_1 has first maximum near x = 1.841
        REQUIRE(SF.J1(1.841) > SF.J1(1.5));
        REQUIRE(SF.J1(1.841) > SF.J1(2.2));
    }
}
