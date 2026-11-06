// Double-double arithmetic (ATOOLS/Math/DDouble.H): the error-free
// transformations the YFS Comix currents rely on for soft and collinear
// propagators. Each test is a case where plain double loses the answer.
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "ATOOLS/Math/DDouble.H"

#include <cmath>

using ATOOLS::DDouble;

TEST_CASE("DDouble keeps what double rounds away", "[ATOOLS][DDouble]") {
  // 1 + 1e-20 is 1 in double; the double-double carries the 1e-20
  const DDouble a(DDouble(1.0) + DDouble(1e-20));
  CHECK(a.High() == 1.0);
  CHECK(a.Low() == 1e-20);
  CHECK((double)(a - DDouble(1.0)) == 1e-20);
  CHECK(a > DDouble(1.0));
  CHECK(a != DDouble(1.0));
}

TEST_CASE("DDouble product is exact", "[ATOOLS][DDouble]") {
  // (1+2^-30)^2 = 1 + 2^-29 + 2^-60: the last term is below double's
  // resolution at 1 and must end up in the low word
  const double e(std::ldexp(1.0, -30));
  const DDouble p(DDouble(1.0 + e) * DDouble(1.0 + e));
  CHECK(p.High() == 1.0 + 2.0*e);
  CHECK(p.Low() == e*e);
}

TEST_CASE("DDouble division and square root", "[ATOOLS][DDouble]") {
  const DDouble third(DDouble(1.0) / DDouble(3.0));
  CHECK(std::abs((double)(third*DDouble(3.0) - DDouble(1.0))) < 1e-30);
  const DDouble r(sqrt(DDouble(2.0)));
  CHECK(std::abs((double)(r*r - DDouble(2.0))) < 1e-30);
  CHECK((double)sqrt(DDouble(0.0)) == 0.0);
  CHECK(isnan(sqrt(DDouble(-1.0))));
}

TEST_CASE("DDouble sign, abs and comparisons", "[ATOOLS][DDouble]") {
  const DDouble x(-2.0, -1e-20);
  CHECK(x < DDouble(-2.0));
  CHECK(abs(x) == DDouble(2.0, 1e-20));
  CHECK(-x == DDouble(2.0, 1e-20));
  CHECK(DDouble(3) == DDouble(3.0));
  CHECK(isfinite(DDouble(1.0)));
}
