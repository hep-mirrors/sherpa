#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Selectors/Selector.H"

#include <memory>

TEST_CASE("FastjetVeto accepts ee points below the jet multiplicity threshold",
          "[PHASIC::Fastjet_Veto]")
{
  PHASIC_Test::SeedGlobals();
  // Only the bunch flavour's lepton classification is needed. Keep it local
  // rather than registering or changing a shared particle-table entry.
  ATOOLS::Particle_Info electron;
  electron.m_kfc = kf_e;
  PHASIC_Test::Bunches bunches(ATOOLS::Flavour(electron),
                               ATOOLS::Flavour(electron, true));

  PHASIC::Selector_Key key;
  bool expected;
  SECTION("zero jets pass a veto on one or more jets") {
    key.m_settings = ATOOLS::Scoped_Settings("FastjetVeto: {N: 1}");
    expected = true;
  }
  SECTION("zero jets fail a veto whose multiplicity threshold is zero") {
    key.m_settings = ATOOLS::Scoped_Settings("FastjetVeto: {N: 0}");
    expected = false;
  }

  std::unique_ptr<PHASIC::Selector_Base> veto(
    PHASIC::Selector_Getter::GetObject("FastjetVeto", key));
  REQUIRE(veto != nullptr);
  ATOOLS::Selector_List particles;
  CHECK(veto->Trigger(particles) == expected);
}
