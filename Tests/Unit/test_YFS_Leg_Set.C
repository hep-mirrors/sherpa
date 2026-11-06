// Leg-set bit operations (YFS/Main/Leg_Set.H): the helicity and leg-subset
// encoding shared by ACRAIC, the NLO reals and the Comix propagator shifts.
// Relabelling helicities between Comix and ACRAIC goes through these.
#include <catch2/catch_test_macros.hpp>

#include "YFS/Main/Leg_Set.H"

#include <set>

using namespace YFS;

TEST_CASE("Leg sets: basic sets", "[YFS][Leg_Set]") {
  CHECK(LegBit(0) == 1u);
  CHECK(LegBit(3) == 8u);
  CHECK(NumLegSubsets(4) == 16u);
  CHECK(AllLegs(4) == 15u);
  CHECK(HasLeg(size_t(10), 1));
  CHECK_FALSE(HasLeg(size_t(10), 2));
  CHECK(LegBitOf(size_t(10), 3) == 1u);
  CHECK(AddLeg(size_t(1), 2) == 5u);
  CHECK(LegUnion(size_t(3), size_t(12)) == 15u);
  CHECK(LegCount(size_t(0b101101)) == 4u);
  CHECK(LowestLeg(size_t(0b101000)) == 3u);
  CHECK(WithoutLowestLeg(size_t(0b101000)) == size_t(0b100000));
  CHECK(RemoveLegs(size_t(15), size_t(5)) == 10u);
}

TEST_CASE("Leg sets: flips are involutions", "[YFS][Leg_Set]") {
  for (size_t cfg(0); cfg < NumLegSubsets(5); ++cfg)
    for (size_t flip(0); flip < NumLegSubsets(5); ++flip)
      CHECK(FlipLegs(FlipLegs(cfg, flip), flip) == cfg);
}

TEST_CASE("Leg sets: Comix packing round trip", "[YFS][Leg_Set]") {
  // Comix stores n fermion helicities first and the photons above them;
  // LowLegs and the shift by n undo JoinLegs
  const size_t n(4), ng(3);
  for (size_t ferm(0); ferm < NumLegSubsets(n); ++ferm)
    for (size_t phot(0); phot < NumLegSubsets(ng); ++phot) {
      const size_t packed(JoinLegs(ferm, phot, n));
      CHECK(LowLegs(packed, n) == ferm);
      CHECK((packed >> n) == phot);
      CHECK(packed < NumLegSubsets(n + ng));
    }
}

TEST_CASE("Leg sets: helicity bits", "[YFS][Leg_Set]") {
  CHECK(HelicityBit(+1) == 0u);
  CHECK(HelicityBit(-1) == 1u);
  const size_t cfg(AddLeg(AddLeg(size_t(0), 1), 3));
  CHECK(LegHelicity(cfg, 0) == +1);
  CHECK(LegHelicity(cfg, 1) == -1);
  CHECK(LegHelicity(cfg, 3) == -1);
}

TEST_CASE("Leg sets: subset enumeration visits every non-empty subset once",
          "[YFS][Leg_Set]") {
  const size_t of(0b10110);
  std::set<size_t> seen;
  for (size_t s(of); s; s = NextSubset(s, of)) {
    CHECK(RemoveLegs(s, of) == 0u);
    CHECK(seen.insert(s).second);
  }
  CHECK(seen.size() == NumLegSubsets(LegCount(of)) - 1);
}
