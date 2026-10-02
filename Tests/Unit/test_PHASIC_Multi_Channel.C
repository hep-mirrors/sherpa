#include <catch2/catch_all.hpp>

#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Channels/Single_Channel.H"

namespace {

class Fixed_Weight_Channel : public PHASIC::Single_Channel {
public:
  double next_weight = 0.;
  void GenerateWeight(const int&) override { SetWeight(next_weight); }
};

} // namespace

TEST_CASE("Multi_Channel preserves a single incoming channel's zero weight",
          "[PHASIC::Multi_Channel]")
{
  PHASIC::Multi_Channel channels("single incoming channel");
  auto* channel = new Fixed_Weight_Channel();
  channels.Add(channel); // Multi_Channel owns the channel.

  SECTION("a vetoed first point has zero weight, not the constructor weight") {
    channels.GenerateWeight();
    CHECK(channels.Weight() == Catch::Approx(0.).margin(1e-15));
  }
  SECTION("a vetoed point does not reuse the previous point's nonzero weight") {
    channel->next_weight = 2.5;
    channels.GenerateWeight();
    REQUIRE(channels.Weight() == Catch::Approx(channel->next_weight));

    channel->next_weight = 0.;
    channels.GenerateWeight();
    CHECK(channels.Weight() == Catch::Approx(0.).margin(1e-15));

    channel->next_weight = 4.;
    channels.GenerateWeight();
    CHECK(channels.Weight() == Catch::Approx(channel->next_weight));
  }
}
