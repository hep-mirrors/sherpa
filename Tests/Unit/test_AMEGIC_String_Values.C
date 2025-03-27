#include <catch2/catch_all.hpp>
#include "AMEGIC++/String/Values.H"

TEST_CASE("Value construction", "[AMEGIC::String::Values]") {
  AMEGIC::Values v;
  CHECK(true);
}

TEST_CASE("Value calculation", "[AMEGIC::String::Values]") {
  AMEGIC::Values v;
  v.Calculate();
  CHECK(true);
}

TEST_CASE("Setting coupling flavours", "[AMEGIC::String::Values]") {
  AMEGIC::Values v;
  std::vector<Complex> cvec{Complex(0.,0.)};
  v.SetCouplFlav(cvec);
  CHECK(true);
}

TEST_CASE("Retrieving number of coupling flavours", "[AMEGIC::String::Values]") {
  AMEGIC::Values v;
  CHECK(!v.NumberOfCouplings());
}

TEST_CASE("Value evaluation", "[AMEGIC::String::Values]") {
  AMEGIC::Values v;
  CHECK(v.Evaluate(0,0) == Complex(0.,0.));
}

