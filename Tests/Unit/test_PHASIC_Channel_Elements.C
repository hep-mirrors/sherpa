#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Channels/Channel_Elements.H"

TEST_CASE("Channel_Elements exponential inversion returns the density and random coordinate",
          "[PHASIC::Channel_Elements]")
{
  PHASIC_Test::SeedGlobals();
  const double lower = 10., upper = 30.;
  double masses[] = {1., 1.}; // Unused by this invariant-mass map.
  for (const double slope : {.02, .1}) {
    for (const double random : {0., .2, .8, 1.}) {
      CAPTURE(slope, random);
      const double value = PHASIC::CE.ExponentialMomenta(
        slope, lower, upper, masses, random);
      double recovered = -1.;
      const double density = PHASIC::CE.ExponentialWeight(
        slope, lower, upper, masses, value, recovered);
      // Every propagator Weight helper supplies a density to channels, which
      // take its reciprocal. Differentiate the independently normalized CDF.
      const double expected = slope * std::exp(-slope * (value - lower))
                              / (1. - std::exp(-slope * (upper - lower)));
      CHECK_THAT(density, Catch::Matchers::WithinRel(expected, 1.e-12));
      CHECK_THAT(recovered, Catch::Matchers::WithinAbs(random, 1.e-12));
    }
  }
}
