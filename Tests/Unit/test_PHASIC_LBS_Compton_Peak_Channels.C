#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Channels/LBS_Compton_Peak_Channels.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <memory>

TEST_CASE("LBS_Compton_Peak channels invert their generated invariant-mass map",
          "[PHASIC::LBS_Compton_Peak_Channels]")
{
  PHASIC_Test::SeedGlobals();
  const double exponent = .5;
  for (const double pole : {.005, .4, .9}) {
    for (const int shape : {0, 1, 2, 3}) {
      CAPTURE(pole, shape);
      ATOOLS::Integration_Info info;
      std::unique_ptr<PHASIC::Single_Channel> channel;
      if (shape == 0)
        channel.reset(new PHASIC::LBS_Compton_Peak_Uniform(
          exponent, pole, "BEAM", &info, 3));
      else if (shape == 1)
        channel.reset(new PHASIC::LBS_Compton_Peak_Central(
          exponent, pole, "BEAM", &info, 3));
      else if (shape == 2)
        channel.reset(new PHASIC::LBS_Compton_Peak_Forward(
          exponent, pole, .5, "BEAM", &info, 3));
      else
        channel.reset(new PHASIC::LBS_Compton_Peak_Backward(
          exponent, pole, .5, "BEAM", &info, 3));
      ATOOLS::Info_Key sp, x, y, recovered;
      sp.Assign("BEAM::s'", 5, 0, &info);
      x.Assign("BEAM::x", 6, 0, &info);
      y.Assign("BEAM::y", 3, 0, &info);
      recovered.Assign("LBS_Compton_Peak_" + ATOOLS::ToString(exponent) +
                       "_" + ATOOLS::ToString(pole), 1, 0, &info);
      for (const double random_mass : {.2, .8}) {
        info.ResetAll();
        sp[0] = 1.; sp[1] = 80.; sp[2] = 100.; sp[4] = 0.;
        x[0] = x[1] = std::log(.001); x[2] = x[3] = 0.;
        y[0] = -10.; y[1] = 10.;
        const double random[] = {random_mass, .3};
        channel->GeneratePoint(random);
        REQUIRE(sp[3] >= sp[0]);
        REQUIRE(sp[3] <= sp[1]);
        channel->GenerateWeight();
        CHECK_THAT(recovered[0], Catch::Matchers::WithinAbs(random_mass, 1.e-12));
        CHECK(channel->Weight() > 0.);
      }
    }
  }
}
