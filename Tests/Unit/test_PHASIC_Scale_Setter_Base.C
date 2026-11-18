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

// sets all scales to hH_T2 of the (sub)event
class HadronicHT : public PHASIC::Scale_Setter_Base {
public:
  HadronicHT(const PHASIC::Scale_Setter_Arguments& args, ATOOLS::NLO_subevtlist* subs)
    : Scale_Setter_Base(args) { p_subs = subs; }
private:
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override
  {
    const double ht = hHT();
    m_scale.assign(ATOOLS::stp::size, ht * ht);
    return m_scale[ATOOLS::stp::fac];
  }
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

TEST_CASE("Flavour-dependent scales use the flavours of each NLO subevent",
          "[PHASIC::Scale_Setter_Base]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle electron(kf_e), up(kf_u), photon(kf_photon);
  up.Info().m_strong = 3;
  const ATOOLS::Flavour e(kf_e), u(kf_u), a(kf_photon);
  // e- e+ -> u ubar gamma at NLO EW in Sherpa's order, the colourless photon
  // before the quarks; the born subevent of the photon emission lacks the
  // photon, so its quarks sit one index lower than in the real process
  PHASIC_Test::Process real({e, e.Bar(), a, u, u.Bar()}, 2);
  PHASIC_Test::Process born({e, e.Bar(), u, u.Bar()}, 2);
  const ATOOLS::Flavour_Vector realfl = real.Flavours(), bornfl = born.Flavours();
  const ATOOLS::Vec4D_Vector realp = {ATOOLS::Vec4D(50., 0., 0., 50.),
                                     ATOOLS::Vec4D(50., 0., 0., -50.),
                                     ATOOLS::Vec4D(10., 0., 6., 8.),
                                     ATOOLS::Vec4D(13., 5., 0., 12.),
                                     ATOOLS::Vec4D(25., -7., 0., -24.)};
  const ATOOLS::Vec4D_Vector bornp = {ATOOLS::Vec4D(50., 0., 0., 50.),
                                     ATOOLS::Vec4D(50., 0., 0., -50.),
                                     ATOOLS::Vec4D(5., 3., 0., 4.),
                                     ATOOLS::Vec4D(5., -4., 0., -3.)};
  ATOOLS::NLO_subevt dipole(4, nullptr, &bornfl.front(), &bornp.front());
  ATOOLS::NLO_subevt realsub(5, nullptr, &realfl.front(), &realp.front());
  dipole.p_proc = static_cast<PHASIC::Process_Base*>(&born);
  realsub.p_proc = static_cast<PHASIC::Process_Base*>(&real);
  dipole.m_trig = realsub.m_trig = true;
  ATOOLS::NLO_subevtlist subs;
  subs.push_back(&dipole);
  subs.push_back(&realsub);

  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = &real;
  HadronicHT setter(args, &subs);
  setter.CalculateScale(realp);
  // hH_T2 sums both quarks, not the photon
  CHECK_THAT(dipole.m_mu2[ATOOLS::stp::fac], Catch::Matchers::WithinRel(49., 1.e-12));
  CHECK_THAT(realsub.m_mu2[ATOOLS::stp::fac], Catch::Matchers::WithinRel(144., 1.e-12));
}
