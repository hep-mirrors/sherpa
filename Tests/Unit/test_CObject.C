#include <catch2/catch_test_macros.hpp>
#include "METOOLS/Currents/C_Spinor.H"

TEST_CASE("Spinor construction", "[spinor]") {
  ATOOLS::Vec4<double> v0(0, 0, 0, 0);
  METOOLS::CSpinor<double> s0(1, 1, 1, v0);
  REQUIRE(s0.IsZero());
}

// ########################################
// for Catch2 v2 use the below code instead

/*
#include <catch2/catch.hpp>
#include "METOOLS/Currents/C_Spinor.H"

TEST_CASE("Spinor construction", "[spinor]") {
  ATOOLS::Vec4<double> v0(0, 0, 0, 0);
  METOOLS::CSpinor<double> s0(1, 1, 1, v0);
  REQUIRE(s0.IsZero());
}
*/

// and run the following bash line here:
/*
echo "#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>" > test_main.C
*/
