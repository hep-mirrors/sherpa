#ifndef TESTS_UNIT_PHASIC_MERGING_SCALE_TEST_UTILS_H
#define TESTS_UNIT_PHASIC_MERGING_SCALE_TEST_UTILS_H

#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Scale_Test_Utils.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/Core_Scale_Setter.H"
#include "PDF/Main/Shower_Base.H"

#include <memory>
#include <stdexcept>

namespace PHASIC_Test {

// Constructor-only tests need a shower handle, not physical clustering or an
// initialized generator. Never call Calculate/Initialize with this stub.
class ConstructorShower : public PDF::Shower_Base {
public:
  ConstructorShower() : PDF::Shower_Base("PHASIC constructor fixture") {}
  int PerformShowers() override { return 0; }
  int PerformDecayShowers() override { return 0; }
  bool ExtractPartons(ATOOLS::Blob_List*) override { return false; }
  void CleanUp() override {}
  PDF::Cluster_Definitions_Base* GetClusterDefinitions() override { return nullptr; }
  bool PrepareShower(ATOOLS::Cluster_Amplitude*, const bool&) override { return false; }
};

inline std::unique_ptr<PHASIC::Scale_Setter_Base> ConstructMergingScale(
    Process& process, const std::string& tag)
{
  PHASIC::Scale_Setter_Arguments args(nullptr, tag);
  args.p_proc = &process;
  return std::unique_ptr<PHASIC::Scale_Setter_Base>(
    PHASIC::Scale_Setter_Base::Scale_Getter_Function::GetObject(tag, args));
}

inline void CheckMalformedScaleHeader(const std::string& name)
{
  ScaleEnvironment environment;
  ConstructorShower shower;
  Process process(ATOOLS::Flavour_Vector(4, ATOOLS::Flavour(kf_none)), 2);
  process.SetShower(&shower);
  REQUIRE(ConstructMergingScale(process, name) != nullptr);
  REQUIRE(ConstructMergingScale(process, name + "[M=2]") != nullptr);
  REQUIRE_THROWS_AS(ConstructMergingScale(process, name + "[M=2"), ATOOLS::fatal_error);
}

class ConstructorCore : public PHASIC::Core_Scale_Setter {
  int& m_live;
public:
  ConstructorCore(const PHASIC::Core_Scale_Arguments& args, int& live)
    : PHASIC::Core_Scale_Setter(args), m_live(live) { ++m_live; }
  ~ConstructorCore() override { --m_live; }
  PDF::Cluster_Param Calculate(ATOOLS::Cluster_Amplitude*) override
  { throw std::logic_error("constructor-only core fixture"); }
};

// Scoped registration avoids leaving a test getter or its counter in the
// global factory after the objects under test have been destroyed.
class ConstructorCoreFactory : public PHASIC::Core_Scale_Getter {
  int& m_live;
  PHASIC::Core_Scale_Setter* operator()(const PHASIC::Core_Scale_Arguments& args) const override
  { return new ConstructorCore(args, m_live); }
public:
  explicit ConstructorCoreFactory(int& live)
    : PHASIC::Core_Scale_Getter("PHASIC_InlineUnorderedCore"), m_live(live) {}
};

inline void CheckInlineUnorderedScale(const std::string& name)
{
  ScaleEnvironment environment;
  ConstructorShower shower;
  Process process(ATOOLS::Flavour_Vector(4, ATOOLS::Flavour(kf_none)), 2);
  process.SetShower(&shower);
  int live = 0;
  ConstructorCoreFactory factory(live);
  {
    auto setter = ConstructMergingScale(process, name + "[U=PHASIC_InlineUnorderedCore]");
    REQUIRE(setter != nullptr);
    CHECK(live == 1);
  }
  CHECK(live == 0);
}

} // namespace PHASIC_Test

#endif
