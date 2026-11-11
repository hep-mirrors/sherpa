#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Channels/FSR_Channel.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Multi_Channel.H"

#include <memory>

TEST_CASE("Decay2Channel uses the full parent/daughter flavour layout",
          "[PHASIC::Decay2Channel]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle parent(kf_h0), first(kf_e), second(kf_mu);
  parent.SetMass(10.);
  first.SetMass(2.);
  second.SetMass(3.);
  const ATOOLS::Flavour_Vector flavours = {ATOOLS::Flavour(kf_h0),
                                           ATOOLS::Flavour(kf_e),
                                           ATOOLS::Flavour(kf_mu)};
  PHASIC_Test::Process process(flavours, 1);
  PHASIC::Multi_Channel channels("decay-test");
  SECTION("direct construction reads both outgoing masses") {
    channels.Add(new PHASIC::Decay2Channel(1, 2, flavours.data()));
  }
  SECTION("the registered factory passes parent then daughters") {
    std::unique_ptr<PHASIC::Channel_Generator> generator(
      PHASIC::Channel_Generator::Getter_Function::GetObject(
        "Decay2", PHASIC::Channel_Generator_Key("Decay2", &process, &channels)));
    REQUIRE(generator != nullptr);
    generator->GenerateChannels();
  }
  REQUIRE(channels.Number() == 1);
  ATOOLS::Vec4D momenta[3] = {ATOOLS::Vec4D(10., 0., 0., 0.)};
  double random[] = {.25, .7};
  auto* channel = channels.Channel(0);
  channel->GeneratePoint(momenta, nullptr, random);
  const auto residual = momenta[0] - momenta[1] - momenta[2];
  for (size_t i = 0; i < 4; ++i)
    CHECK_THAT(residual[i], Catch::Matchers::WithinAbs(0., 1.e-12));
  CHECK_THAT(momenta[1].Abs2(), Catch::Matchers::WithinAbs(4., 1.e-12));
  CHECK_THAT(momenta[2].Abs2(), Catch::Matchers::WithinAbs(9., 1.e-12));
  CHECK_THAT(momenta[1][0], Catch::Matchers::WithinAbs((100. + 4. - 9.) / 20., 1.e-12));
  channel->GenerateWeight(momenta, nullptr);
  const double phase_volume = std::sqrt(std::pow(100. - 4. - 9., 2) - 4. * 4. * 9.)
                              / 100. / (8. * M_PI);
  CHECK_THAT(channel->Weight(), Catch::Matchers::WithinRel(phase_volume, 1.e-12));
}
