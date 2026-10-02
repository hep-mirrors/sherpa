#include <catch2/catch_all.hpp>

#include "PHASIC++/Channels/Vegas.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"

#include <vector>

namespace {

void SeedVegasGlobals()
{
  if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
  if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
  static bool settings_seeded = false;
  if (!settings_seeded) {
    ATOOLS::Settings::InitializeMainSettings(std::string{"VEGAS_MODE: 2"});
    settings_seeded = true;
  }
}

} // namespace

TEST_CASE("Vegas refinement preserves the map and its Jacobian",
          "[PHASIC::Vegas]")
{
  SeedVegasGlobals();
  int dimension = 1;
  bool preset = false;
  SECTION("splitting a nonuniform one-dimensional grid preserves sampling") {
    preset = true;
  }
  SECTION("a dimension count larger than the refined bin count is safe") {
    dimension = 21; // Default adaptive mode refines 10 bins into 20.
  }

  PHASIC::Vegas vegas(dimension, 100, "refinement_regression");
  PHASIC::Vegas::Power mapping(1.);
  if (preset) vegas.Preset(&mapping, 0);
  const int initial_bins = vegas.GetNBins();
  if (dimension > 1) REQUIRE(dimension > 2 * initial_bins);

  const std::vector<double> samples{0.025, 0.37, 0.73, 0.975};
  std::vector<std::vector<double>> before;
  std::vector<double> weights;
  for (const double sample : samples) {
    std::vector<double> random(dimension, sample);
    const double* point = vegas.GeneratePoint(random.data());
    before.emplace_back(point, point + dimension);
    weights.push_back(vegas.GenerateWeight(point));
  }

  vegas.Refine();
  CHECK(vegas.GetNBins() == 2 * initial_bins);
  for (size_t sample = 0; sample < samples.size(); ++sample) {
    std::vector<double> random(dimension, samples[sample]);
    const double* point = vegas.GeneratePoint(random.data());
    for (int i = 0; i < dimension; ++i)
      CHECK(point[i] == Catch::Approx(before[sample][i]).epsilon(1e-12));
    CHECK(vegas.GenerateWeight(point)
          == Catch::Approx(weights[sample]).epsilon(1e-12));
  }
}
