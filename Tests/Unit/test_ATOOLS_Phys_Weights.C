#include <catch2/catch_all.hpp>

#include "ATOOLS/Phys/Weights.H"

#include <cmath>
#include <limits>
#include <string>

using ATOOLS::Weights;
using ATOOLS::Weights_Map;

namespace {

  // Build the relative weights map of a single event, with the matrix-element
  // variations `me` stored in both "Main" (ME-only) and "All" (ME & PS), and
  // a variation-independent shower weight `sud` multiplied in, as done by the
  // shower when its own reweighting is disabled.
  Weights_Map EventWeights(double base, double var1, double var2, double sud,
                           bool sud_in_main)
  {
    Weights me(1.0);
    me["v1"] = var1;
    me["v2"] = var2;
    Weights_Map ev(base);
    ev["Main"] = me;
    ev["All"] = me;
    ev["Sudakov"] = Weights(sud);
    ev["All"] *= sud;
    if (sud_in_main)
      ev["Main"] *= sud;
    return ev;
  }

} // namespace

// Regression test for #667: the nominal shower weight (e.g. the MC@NLO
// accept/reject weight) differs from event to event. Unless it is also
// carried by "Main", summed ME-only variations are not reweighted consistently
// with the nominal, and deviate from the ME & PS ones even with shower
// reweighting disabled.
TEST_CASE("ME-only and ME & PS variation sums agree without shower "
          "reweighting",
          "[ATOOLS::Weights_Map]") {
  const bool sud_in_main = GENERATE(true, false);

  Weights_Map sum(0.0);
  sum += EventWeights(2.0, 1.5, 0.5, 3.0, sud_in_main);
  sum += EventWeights(1.0, 0.8, 1.2, -0.5, sud_in_main);
  sum.MakeAbsolute();

  const Weights& main = sum.at("Main");
  const Weights& all = sum.at("All");
  REQUIRE(main.Size() == 3);
  REQUIRE(all.Size() == 3);

  // ME & PS: 2*{1, 1.5, 0.5}*3 + 1*{1, 0.8, 1.2}*(-0.5)
  CHECK_THAT(all[0], Catch::Matchers::WithinRel(5.5));
  CHECK_THAT(all[1], Catch::Matchers::WithinRel(8.6));
  CHECK_THAT(all[2], Catch::Matchers::WithinRel(2.4));

  if (sud_in_main) {
    for (size_t i {0}; i < 3; ++i)
      CHECK_THAT(main[i], Catch::Matchers::WithinRel(all[i]));
  } else {
    // Rescaling the ME-only sums to the nominal does not recover ME & PS.
    for (size_t i {1}; i < 3; ++i)
      CHECK_FALSE(main[i] * all[0] / main[0] ==
                  Catch::Approx(all[i]).epsilon(1e-3));
  }
}

// An MC@NLO H-event in which R and DS cancel exactly in the nominal and in all
// variations but one, which keeps the smallest subnormal as round-off (as
// found in a photoproduction NLO+PS run, about one event in 5e5). With all
// nominals zero, MakeRelative normalises by that subnormal, and must not turn
// the exact zeros into NaN (0 * (1/denorm_min) = 0 * inf).
TEST_CASE("MakeRelative stays finite when the only non-zero entry is "
          "subnormal",
          "[ATOOLS::Weights_Map]") {
  const double tiny {-std::numeric_limits<double>::denorm_min()};

  Weights me(1.0);
  me["v1"] = 1.0;
  me["v2"] = 1.0;
  me["v3"] = 1.0;
  Weights_Map ev(1.0);
  ev["Main"] = me;
  ev["All"] = me;
  ev.MakeAbsolute();
  for (const std::string key : {"Main", "All"}) {
    Weights& w = ev[key];
    for (size_t i {0}; i < w.Size(); ++i)
      w[i] = 0.0;
    w[3] = tiny;
  }
  ev.MakeRelative();

  CHECK(std::isfinite(ev.BaseWeight()));
  CHECK(ev.Nominal() == 0.0);
  // nothing but round-off: a zero weight, which the event generation discards
  CHECK(ev.IsZero());
  ev.MakeAbsolute();
  for (const std::string key : {"Main", "All"}) {
    const Weights& w = ev.at(key);
    for (size_t i {0}; i < w.Size(); ++i) {
      INFO(key << "[" << i << "] = " << w[i]);
      CHECK(std::isfinite(w[i]));
      CHECK(std::abs(w[i]) <= std::abs(tiny));
    }
  }
}
