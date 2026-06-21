#include <catch2/catch_all.hpp>

#include "AHADIC++/Tools/Proto_Particle.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

// Regression for commit f12981f0d: the Proto_Particle Flavour& constructor (the
// default-args, most-used ctor) used to leave both m_gen and m_kt2max
// uninitialised. Pin that it now initialises generation to 1 and kt2max to the
// energy squared. The ctor stores the flavour without interrogating it, so no KF
// table or hadpars singleton is needed.
//
// The ATOOLS::Particle& ctor (which the same fix also touched, for m_gen) is not
// covered here: constructing an ATOOLS::Particle pulls in global state that is
// not proportionate to seed in a unit test.
TEST_CASE("Proto_Particle Flavour& ctor initialises generation and kt2max",
          "[AHADIC::Proto_Particle]") {
  if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();

  SECTION("generation is 1 and kt2max is E^2") {
    Proto_Particle p(Flavour(kf_d), Vec4D(2., 0., 0., 1.));
    CHECK(p.Generation() == 1);
    CHECK_THAT(p.KT2_Max(), Catch::Matchers::WithinRel(sqr(2.), 1.e-12));
  }
}
