// for Catch2 v3 change <catch2/catch.hpp> to <catch2/catch_test_macros.hpp>
// and remove test_main.C
#include <catch2/catch.hpp>
#include "METOOLS/Currents/C_Spinor.H"

TEST_CASE("Spinor construction", "[spinor]") {
  ATOOLS::Vec4<double> v0(0, 0, 0, 0);
  METOOLS::CSpinor<double> s0(1, 1, 1, v0);
  REQUIRE(s0.IsZero());
}
