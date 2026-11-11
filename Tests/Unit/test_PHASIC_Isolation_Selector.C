#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Selectors/Selector.H"

#include <memory>

namespace {

class Lifetime_Probe : public PHASIC::Selector_Base {
  size_t& m_destroyed;
public:
  explicit Lifetime_Probe(size_t& destroyed)
    : Selector_Base("IsolationLifetimeProbe"), m_destroyed(destroyed) {}
  ~Lifetime_Probe() override { ++m_destroyed; }
  bool Trigger(ATOOLS::Selector_List&) override { return true; }
  void BuildCuts(PHASIC::Cut_Data*) override {}
};

class Lifetime_Probe_Getter : public PHASIC::Selector_Getter {
  size_t& m_created;
  size_t& m_destroyed;
public:
  Lifetime_Probe_Getter(size_t& created, size_t& destroyed)
    : PHASIC::Selector_Getter("IsolationLifetimeProbe"),
      m_created(created), m_destroyed(destroyed) {}
protected:
  PHASIC::Selector_Base* operator()(const PHASIC::Selector_Key&) const override
  {
    ++m_created;
    return new Lifetime_Probe(m_destroyed);
  }
};

} // namespace

TEST_CASE("Isolation_Selector destroys its owned subselectors",
          "[PHASIC::Isolation_Selector]")
{
  // The null-process fixture has no external legs. Only kf_none is needed
  // by Isolation_Selector's default flavour and Output_ID.
  PHASIC_Test::SeedGlobals();
  size_t created = 0, destroyed = 0;
  Lifetime_Probe_Getter getter(created, destroyed);
  PHASIC::Selector_Key key;

  SECTION("two direct children are each destroyed once") {
    key.m_settings = ATOOLS::Scoped_Settings(R"(
Isolation_Selector:
  Subselectors:
    - [IsolationLifetimeProbe]
    - [IsolationLifetimeProbe]
)");
  }
  SECTION("destruction recurses through nested isolation selectors") {
    key.m_settings = ATOOLS::Scoped_Settings(R"(
Isolation_Selector:
  Subselectors:
    - Isolation_Selector:
        Subselectors:
          - [IsolationLifetimeProbe]
          - [IsolationLifetimeProbe]
)");
  }

  {
    std::unique_ptr<PHASIC::Selector_Base> selector(
      PHASIC::Selector_Getter::GetObject("Isolation_Selector", key));
    REQUIRE(selector != nullptr);
    REQUIRE(created == 2);
    CHECK(destroyed == 0);
  }
  CHECK(destroyed == created);
}
