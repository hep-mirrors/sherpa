#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include <array>
#include <memory>
#include <new>

namespace {
class Scale : public PHASIC::Scale_Setter_Base {
public:
  explicit Scale(const PHASIC::Scale_Setter_Arguments& args) : Scale_Setter_Base(args) {}
private:
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override { return 0.; }
};
}

TEST_CASE("Scale setter arguments initialize fallback multiplicities",
          "[PHASIC::Scale_Setter_Base]")
{
  PHASIC_Test::SeedGlobals();
  using Arguments = PHASIC::Scale_Setter_Arguments;
  alignas(Arguments) std::array<unsigned char, sizeof(Arguments)> storage;
  storage.fill(0x3f); // Do not let allocator zeroing mask unspecified members.
  auto destroy = [](Arguments* args) { args->~Arguments(); };
  std::unique_ptr<Arguments, decltype(destroy)> args(new (storage.data()) Arguments(nullptr), destroy);
  REQUIRE(args->m_nin == 0);
  REQUIRE(args->m_nout == 0);
  SECTION("without a process the default setter has no external legs") {
    Scale setter(*args);
    CHECK(setter.NIn() == 0);
    CHECK(setter.NOut() == 0);
    CHECK(setter.Momenta().empty());
  }
  SECTION("explicit multiplicities are used without a process") {
    args->m_nin = 2;
    args->m_nout = 3;
    Scale setter(*args);
    CHECK(setter.NIn() == 2);
    CHECK(setter.NOut() == 3);
    CHECK(setter.Momenta().size() == 5);
  }
  SECTION("process metadata takes precedence over fallback multiplicities") {
    PHASIC_Test::Process process(ATOOLS::Flavour_Vector(4, ATOOLS::Flavour(kf_none)), 1);
    args->p_proc = &process;
    Scale setter(*args);
    CHECK(setter.NIn() == 1);
    CHECK(setter.NOut() == 3);
    CHECK(setter.Momenta().size() == 4);
  }
}
