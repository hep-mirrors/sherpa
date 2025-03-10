#include <catch2/catch_test_macros.hpp>
#include "AMEGIC++/String/Values.H"

TEST_CASE("AMEGIC Values class", "[values]") {
  AMEGIC::Values v;
  std::vector<Complex> cvec{Complex(0.,0.)};
  v.Calculate();
  v.SetCouplFlav(cvec);
  const bool result = !v.NumberOfCouplings() && v.Evaluate(0,0) == Complex(0.,0.);
  REQUIRE(result);
}

