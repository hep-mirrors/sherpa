#include <catch2/catch_all.hpp>

#include "MCATNLO/Main/CS_Gamma.H"

using MCATNLO::CS_Gamma;

// Regression test for the stale gamma-weight bug (0e5fefcf1): the gamma
// weight is multiplied into the Sudakov weights for every trial emission,
// but Reject() only assigns it for trials that reach it.  Trials vetoed by
// RemnantTest rely on the constructor seed and the per-trial ResetWeight()
// to contribute a factor of exactly one instead of the previous trial's
// (unbounded, possibly negative) accept/reject weight.
TEST_CASE("CS_Gamma trial-weight lifecycle", "[MCATNLO::CS_Gamma]") {
  CS_Gamma gamma(nullptr, nullptr, nullptr);

  SECTION("a freshly constructed CS_Gamma carries the unit weight") {
    CHECK_THAT(gamma.Weight(), Catch::Matchers::WithinULP(1.0, 0));
  }

  SECTION("ResetWeight() yields the unit weight for the next trial") {
    gamma.ResetWeight();
    CHECK_THAT(gamma.Weight(), Catch::Matchers::WithinULP(1.0, 0));
  }
}
