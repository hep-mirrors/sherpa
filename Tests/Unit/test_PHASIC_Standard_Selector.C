#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"

#include <memory>

TEST_CASE("Polar_Angle preserves explicit angle units when a frame flag is given",
          "[PHASIC::Standard_Selector]")
{
  PHASIC_Test::SeedGlobals();
  if (ATOOLS::s_kftable.find(kf_e) == ATOOLS::s_kftable.end())
    ATOOLS::AddParticle(kf_e, 0.000511, 0., 0., -3, 1, true, 1, "e-", "e");
  ATOOLS::Flavour flavour(kf_e);
  PHASIC_Test::Process process({flavour, flavour, flavour}, 2);

  for (const std::string parameters : {
         "1, 2", "80, 100, 0", "80, 100, 0, 0", "1, 2, 1, 0"}) {
    CAPTURE(parameters);
    PHASIC::Selector_Key key;
    key.p_proc = &process;
    key.m_settings = ATOOLS::Scoped_Settings(
      "[Polar_Angle, 11, " + parameters + "]");
    std::unique_ptr<PHASIC::Selector_Base> selector(
      PHASIC::Selector_Getter::GetObject("Polar_Angle", key));
    REQUIRE(selector != nullptr);
    ATOOLS::Selector_List particles(process.Flavours(),
      {ATOOLS::Vec4D(10., 0., 0., 10.),
       ATOOLS::Vec4D(10., 0., 0., -10.),
       ATOOLS::Vec4D(1., 1., 0., 0.)}, 2);
    CHECK(selector->Trigger(particles)); // pi/2 radians = 90 degrees
    particles[2].SetMomentum(
      ATOOLS::Vec4D(1., std::sin(0.3), 0., std::cos(0.3)));
    CHECK_FALSE(selector->Trigger(particles));
  }
}
