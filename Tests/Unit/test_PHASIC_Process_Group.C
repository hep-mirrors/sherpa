#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Process/Process_Group.H"
#include "PHASIC++/Main/Phase_Space_Point.H"

namespace {

class Group : public PHASIC::Process_Group {
public:
  PHASIC::Process_Base* GetProcess(const PHASIC::Process_Info&) const override
  { return nullptr; }
  bool Initialize(PHASIC::Process_Base*) override { return true; }
  PHASIC::Process_Base* Find(const std::string& name)
  { return GetProcess(name); }
private:
  using PHASIC::Process_Group::GetProcess;
};

class Child : public PHASIC_Test::Process {
  size_t& m_destroyed;
public:
  Child(const std::string& name, size_t& destroyed)
    : PHASIC_Test::Process({}, 0), m_destroyed(destroyed)
  { m_name = name; }
  ~Child() override { ++m_destroyed; }
  int weight_mode = -1, phase_space_mode = -1;
  ATOOLS::Variations_Mode variation_mode = ATOOLS::Variations_Mode::all;
  ATOOLS::Weight_Info* OneEvent(int wmode, ATOOLS::Variations_Mode varmode,
                               int mode) override
  {
    weight_mode = wmode;
    variation_mode = varmode;
    phase_space_mode = mode;
    return nullptr;
  }
};

} // namespace

TEST_CASE("Process_Group::Clear releases children and resets selection and lookups",
          "[PHASIC::Process_Group]")
{
  PHASIC_Test::SeedGlobals();
  size_t destroyed = 0;
  {
    Group group;
    auto* first = new Child("first", destroyed);
    auto* second = new Child("second", destroyed);
    group.Add(first);
    group.Add(second);
    REQUIRE(group.Find("first") == first);
    REQUIRE(group.Find("second") == second);
    SECTION("the group has not selected a child") {}
    SECTION("the selected child is deleted before the other child") {
      REQUIRE(group.SetSelected(second));
    }
    group.Clear();
    CHECK(destroyed == 2);
    CHECK(group.Size() == 0);
    CHECK(group.Find("first") == nullptr);
    CHECK(group.Find("second") == nullptr);
    CHECK(group.Selected() == &group);
    group.Clear();
    CHECK(destroyed == 2);
    // A cleared group can be repopulated with the same process name.
    group.Add(new Child("first", destroyed));
    CHECK(group.Size() == 1);
  }
  CHECK(destroyed == 3);
}

TEST_CASE("Process_Group forwards all event-generation modes to its child",
          "[PHASIC::Process_Group]")
{
  PHASIC_Test::SeedGlobals();
  // The uniform fallback selection draws from ATOOLS::ran.
  PHASIC_Test::ScopedRandom random;
  size_t destroyed = 0;
  Group group;
  auto* child = new Child("only", destroyed);
  group.Add(child);
  // Zero total cross section uses the uniform fallback selection branch.
  for (const auto varmode : {ATOOLS::Variations_Mode::nominal_only,
                            ATOOLS::Variations_Mode::all}) {
    CHECK(group.OneEvent(1, varmode, PHASIC::psmode::no_gen_isr) == nullptr);
    CHECK(child->weight_mode == 1);
    CHECK(child->variation_mode == varmode);
    CHECK(child->phase_space_mode == PHASIC::psmode::no_gen_isr);
  }
}
