#include <catch2/catch_all.hpp>

#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Terminator_Objects.H"

#include <cmath>
#include <vector>

using namespace ATOOLS;
using namespace PHASIC;

namespace {

  void Boot() {
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::exh) ATOOLS::exh = new ATOOLS::Terminator_Object_Handler();
    if (!ATOOLS::ran) ATOOLS::ran = new ATOOLS::Random(1234);
  }

  // incoming partons with x1 != x2, i.e. a partonic c.m. frame that moves
  // along the beam axis, as produced by the ISR handler in hadron collisions
  std::vector<Vec4D> Incoming(double x1, double x2) {
    const double E = 6800.;
    return {Vec4D(x1*E, 0., 0., x1*E), Vec4D(x2*E, 0., 0., -x2*E)};
  }

  void CheckMomentumConservation(const std::vector<Vec4D>& p, size_t nin,
                                 const std::vector<double>& masses) {
    Vec4D pin, pout;
    for (size_t i(0); i < nin; ++i) pin += p[i];
    for (size_t i(nin); i < p.size(); ++i) {
      pout += p[i];
      CHECK(p[i].Mass() == Catch::Approx(masses[i]).margin(1e-6*pin[0]));
    }
    for (size_t mu(0); mu < 4; ++mu)
      CHECK(pout[mu] == Catch::Approx(pin[mu]).margin(1e-9*pin[0]));
  }

}

TEST_CASE("Rambo generates the final state in the frame of the incoming momenta",
          "[PHASIC][Rambo]") {
  Boot();

  SECTION("2 -> 2 massless, boosted partonic c.m. frame") {
    std::vector<double> masses{0., 0., 0., 0.};
    Rambo rambo(2, masses);
    for (int n(0); n < 100; ++n) {
      std::vector<Vec4D> p(Incoming(0.3, 0.01));
      p.resize(masses.size());
      rambo.GeneratePoint(&p.front(), nullptr);
      CheckMomentumConservation(p, 2, masses);
    }
  }

  SECTION("2 -> 3 massive, boosted partonic c.m. frame") {
    std::vector<double> masses{0., 0., 173., 173., 125.};
    Rambo rambo(2, masses);
    for (int n(0); n < 100; ++n) {
      std::vector<Vec4D> p(Incoming(0.02, 0.4));
      p.resize(masses.size());
      rambo.GeneratePoint(&p.front(), nullptr);
      CheckMomentumConservation(p, 2, masses);
    }
  }

  SECTION("1 -> 3 decay of a moving particle") {
    std::vector<double> masses{173., 80.4, 0., 0.};
    Rambo rambo(1, masses);
    for (int n(0); n < 100; ++n) {
      std::vector<Vec4D> p{Vec4D(std::sqrt(173.*173.+300.*300.), 100., -200., 200.)};
      p.resize(masses.size());
      rambo.GeneratePoint(&p.front(), nullptr);
      CheckMomentumConservation(p, 1, masses);
    }
  }
}

TEST_CASE("The Rambo weight does not depend on the frame of the momenta",
          "[PHASIC][Rambo]") {
  Boot();
  std::vector<double> masses{0., 0., 173., 173., 125.};
  Rambo rambo(2, masses);
  for (int n(0); n < 100; ++n) {
    std::vector<Vec4D> p(Incoming(0.02, 0.4));
    p.resize(masses.size());
    rambo.GeneratePoint(&p.front(), nullptr);
    rambo.GenerateWeight(&p.front(), nullptr);
    const double wlab(rambo.Weight());
    Poincare cms(p[0]+p[1]);
    for (Vec4D& q : p) cms.Boost(q);
    rambo.GenerateWeight(&p.front(), nullptr);
    CHECK(wlab == Catch::Approx(rambo.Weight()).epsilon(1e-9));
  }
}
