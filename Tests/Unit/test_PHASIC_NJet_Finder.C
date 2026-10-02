#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"

#include <cmath>
#include <memory>

TEST_CASE("NJetFinder applies its rapidity window symmetrically",
          "[PHASIC::NJet_Finder]")
{
  PHASIC_Test::SeedGlobals();
  ATOOLS::Particle_Info gluon;
  gluon.m_kfc = kf_gluon;
  gluon.m_strong = 8;
  ATOOLS::Flavour flavour(gluon);
  PHASIC_Test::Process process({flavour, flavour, flavour}, 2);
  PHASIC::Selector_Key key;
  key.p_proc = &process;
  key.m_settings = ATOOLS::Scoped_Settings(
    "NJetFinder: {N: 1, PTMin: 10, YMax: 1}");
  std::unique_ptr<PHASIC::Selector_Base> selector(
    PHASIC::Selector_Getter::GetObject("NJetFinder", key));
  REQUIRE(selector != nullptr);

  for (double y : {-3., -0.5, 0.5, 3.}) {
    CAPTURE(y);
    const double pt = 20.;
    ATOOLS::Selector_List particles(process.Flavours(),
      {ATOOLS::Vec4D(500., 0., 0., 500.),
       ATOOLS::Vec4D(500., 0., 0., -500.),
       ATOOLS::Vec4D(pt * std::cosh(y), pt, 0., pt * std::sinh(y))}, 2);
    CHECK(selector->Trigger(particles) == (std::abs(y) < 1.));
  }
}
