#include <catch2/catch_all.hpp>

#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <memory>

namespace {

struct Bunch_Guard {
  ATOOLS::Flavour first = ATOOLS::rpa->gen.Bunch(0);
  ATOOLS::Flavour second = ATOOLS::rpa->gen.Bunch(1);
  ~Bunch_Guard()
  {
    ATOOLS::rpa->gen.SetBunch(first, 0);
    ATOOLS::rpa->gen.SetBunch(second, 1);
  }
};

} // namespace

TEST_CASE("FastjetVeto accepts ee points below the jet multiplicity threshold",
          "[PHASIC::Fastjet_Veto]")
{
  if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
  if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
  static bool settings_seeded = false;
  if (!settings_seeded) {
    ATOOLS::Settings::InitializeMainSettings(std::string{"{}"});
    settings_seeded = true;
  }

  // Only the bunch flavour's lepton classification is needed. Keep it local
  // rather than registering or changing a shared particle-table entry.
  ATOOLS::Particle_Info electron;
  electron.m_kfc = kf_e;
  Bunch_Guard guard;
  ATOOLS::rpa->gen.SetBunch(ATOOLS::Flavour(electron), 0);
  ATOOLS::rpa->gen.SetBunch(ATOOLS::Flavour(electron, true), 1);

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
