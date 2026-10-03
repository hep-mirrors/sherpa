#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Channels/BBar_Emission_Generator.H"

#include <array>
#include <memory>
#include <new>

namespace {
void SeedSettings()
{
  PHASIC_Test::SeedGlobals();
  ATOOLS::Settings::GetMainSettings()["DIPOLES"]["AMIN"].SetDefault(
    ATOOLS::Max(ATOOLS::rpa->gen.Accu(), 1.e-8));
}
}

TEST_CASE("BBar_Emission_Generator starts without an active emission or weight",
          "[PHASIC::BBar_Emission_Generator]")
{
  SeedSettings();
  using Generator = PHASIC::BBar_Emission_Generator;
  // Pattern ordinary placement-new storage so allocator zeroing cannot hide
  // missing constructor initialization. Do not inspect any private members.
  alignas(Generator) std::array<unsigned char, sizeof(Generator)> storage;
  storage.fill(0x3f);
  auto destroy = [](Generator* generator) { generator->~Generator(); };
  std::unique_ptr<Generator, decltype(destroy)> generator(
    new (storage.data()) Generator(), destroy);
  REQUIRE_THAT(generator->Weight(), Catch::Matchers::WithinAbs(0., 1.e-14));
  const auto active = generator->Active(nullptr);
  CHECK(active.p_dip == nullptr);
  CHECK(active.p_p == nullptr);
  CHECK(active.m_procs.empty());
  CHECK_THAT(active.m_weight, Catch::Matchers::WithinAbs(0., 1.e-14));
  CHECK_FALSE(generator->GenerateWeight(nullptr));
}
