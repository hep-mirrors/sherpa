#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Channel_Elements.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Single_Channel.H"

#include <memory>

TEST_CASE("VHAAG threshold factory preserves fractional thresholds in shared channels",
          "[PHASIC::VHAAG_Threshold]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::Process process(ATOOLS::Flavour_Vector(6, ATOOLS::Flavour(kf_none)), 2);
  PHASIC::Multi_Channel channels("threshold-test");
  std::unique_ptr<PHASIC::Channel_Generator> generator(
    PHASIC::Channel_Generator::Getter_Function::GetObject(
      "VHAAG_Threshold", PHASIC::Channel_Generator_Key(
        "VHAAG_Threshold[I=2,J=3,T=40.5]", &process, &channels)));
  REQUIRE(generator != nullptr);
  generator->GenerateChannels();
  REQUIRE(channels.Number() == 12);

  double random[] = {.3, .4, .2, .6, .7, .35, .25, .8};
  const double expected = PHASIC::CE.ThresholdMomenta(1., 40.5, 0., 10000., random[5]);
  for (int i = 0; i < channels.Number(); ++i) {
    CAPTURE(i);
    ATOOLS::Vec4D momenta[6] = {ATOOLS::Vec4D(50., 0., 0., 50.),
                               ATOOLS::Vec4D(50., 0., 0., -50.)};
    auto* channel = channels.Channel(i);
    channel->GeneratePoint(momenta, nullptr, random);
    CHECK_THAT((momenta[2] + momenta[3]).Abs2(),
               Catch::Matchers::WithinRel(expected, 1.e-10));
    auto residual = momenta[0] + momenta[1];
    for (int j = 2; j < 6; ++j) residual -= momenta[j];
    for (size_t j = 0; j < 4; ++j)
      CHECK_THAT(residual[j], Catch::Matchers::WithinAbs(0., 1.e-10));
    channel->GenerateWeight(momenta, nullptr);
    CHECK(std::isfinite(channel->Weight()));
    CHECK(channel->Weight() > 0.);
  }
}
