#include <catch2/catch_all.hpp>

#include "METOOLS/Explicit/Dipole_Terms.H"
#include "METOOLS/Explicit/Dipole_Info.H"
#include "ATOOLS/Math/Vector.H"

#include <cmath>

using ATOOLS::Vec4D;
using Catch::Matchers::WithinRel;
using namespace METOOLS;

namespace {

  Dipole_Info Info(double alpha) {
    Dipole_Info info;
    info.SetMassive(0);
    for (int t(0); t < 4; ++t) info.SetAMax(t, alpha);
    return info;
  }

  // massless emitter-spectator pair; initial-state legs carry negative
  // energy, so I_Args::type is 0 FF, 1 IF, 2 FI, 3 II
  I_Args Pair(bool is_ij, bool is_k) {
    const Vec4D pij((is_ij ? -1. : 1.) * Vec4D(10., 0., 0., 10.));
    const Vec4D pk((is_k ? -1. : 1.) * Vec4D(10., 0., 0., -10.));
    return I_Args(pij, pk, 0., 0.);
  }

}

TEST_CASE("Only final-final dipoles carry alpha terms in Comix's I-operator",
          "[METOOLS::Dipole_Terms]") {
  const double a = 0.1;
  const Dipole_Info info(Info(a));
  const double la(std::log(a)), ac(a - 1. - std::log(a));

  SECTION("final-final: Nagy's -ln^2(alpha)/2 and collinear (alpha-1-ln alpha) terms") {
    const I_Args ff(Pair(false, false));
    REQUIRE(ff.type == 0);
    CHECK_THAT(FFAE(ff, &info), WithinRel(-0.5 * la * la, 1e-12));
    CHECK_THAT(FFACQQ(ff, &info), WithinRel(1.5 * ac, 1e-12));
    CHECK_THAT(FFACGQ(ff, &info, 0.), WithinRel(-2. / 3. * ac, 1e-12));
    CHECK_THAT(FFACGG(ff, &info), WithinRel(11. / 6. * ac, 1e-12));
  }

  SECTION("an initial-state leg: none, the alpha dependence comes with the KP terms") {
    for (const auto &legs : {std::make_pair(true, false), std::make_pair(false, true),
                             std::make_pair(true, true)}) {
      const I_Args ia(Pair(legs.first, legs.second));
      CAPTURE(ia.type);
      CHECK(FFAE(ia, &info) == 0.);
      CHECK(FFACQQ(ia, &info) == 0.);
      CHECK(FFACGQ(ia, &info, 0.) == 0.);
      CHECK(FFACGG(ia, &info) == 0.);
    }
  }
}
