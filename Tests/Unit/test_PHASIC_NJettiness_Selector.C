#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"

#include <memory>

TEST_CASE("NJettiness rejects unusable clustering algorithms at setup",
          "[PHASIC::NJettiness_Selector]")
{
  // genkt needs an exponent, which NJettiness cannot configure.
  const auto algorithm = GENERATE(as<std::string>{}, "nonexistent", "genkt");
  CAPTURE(algorithm);
  PHASIC_Test::SeedGlobals();
  ATOOLS::Particle_Info gluon;
  gluon.m_kfc = kf_gluon;
  gluon.m_strong = 8;
  ATOOLS::Flavour flavour(gluon);
  PHASIC_Test::Jet_Group jet_group(flavour);
  PHASIC_Test::Bunches bunches(flavour, flavour);
  PHASIC_Test::Process process({flavour, flavour, flavour, flavour}, 2);
  PHASIC::Selector_Key key;
  key.p_proc = &process;
  key.m_settings = ATOOLS::Scoped_Settings(
    "NJettiness: {N: 1, Algorithm: " + algorithm + "}");
  REQUIRE_THROWS_AS(
    std::unique_ptr<PHASIC::Selector_Base>(
      PHASIC::Selector_Getter::GetObject("NJettiness", key)),
    ATOOLS::fatal_error);
}

TEST_CASE("NJettiness uses the requested kt or antikt hard directions",
          "[PHASIC::NJettiness_Selector]")
{
  PHASIC_Test::SeedGlobals();
  ATOOLS::Particle_Info gluon;
  gluon.m_kfc = kf_gluon;
  gluon.m_strong = 8;
  ATOOLS::Flavour flavour(gluon);
  PHASIC_Test::Jet_Group jet_group(flavour);
  PHASIC_Test::Bunches bunches(flavour, flavour);
  PHASIC_Test::Process process(
    {flavour, flavour, flavour, flavour, flavour, flavour}, 2);

  // kt merges the soft parton at phi=0.7 with its nearest neighbour at
  // phi=1.2. antikt instead merges it with the harder parton at phi=0.
  const ATOOLS::Vec4D a(50., 50., 0., 0.);
  const ATOOLS::Vec4D b(10., 10. * std::cos(1.2), 10. * std::sin(1.2), 0.);
  const ATOOLS::Vec4D soft(1., std::cos(0.7), std::sin(0.7), 0.);
  const ATOOLS::Vec4D sum = a + b + soft;
  const ATOOLS::Vec4D recoil(sum.PPerp(), -sum[1], -sum[2], 0.);
  const double beam_e = (sum.E() + recoil.E()) / 2.;
  ATOOLS::Selector_List particles(process.Flavours(),
    {ATOOLS::Vec4D(beam_e, 0., 0., beam_e),
     ATOOLS::Vec4D(beam_e, 0., 0., -beam_e), a, b, soft, recoil}, 2);

  // For jet2 normalization, a two-parton axis contributes M^2/E.
  // In antikt the soft parton still minimizes against the unmerged b axis.
  const double tau_kt = 2. * 10. / 11. * (1. - std::cos(0.5));
  const double tau_antikt = 50. / 51. * (1. - std::cos(0.7))
                          + (1. - std::cos(0.5));
  REQUIRE(tau_kt < tau_antikt);
  const double cut = (tau_kt + tau_antikt) / 2.;
  for (const std::string algorithm : {"kt", "antikt"}) {
    CAPTURE(algorithm);
    PHASIC::Selector_Key key;
    key.p_proc = &process;
    key.m_settings = ATOOLS::Scoped_Settings(
      "NJettiness: {N: 3, TauN: " + std::to_string(cut)
      + ", Algorithm: " + algorithm + "}");
    std::unique_ptr<PHASIC::Selector_Base> selector(
      PHASIC::Selector_Getter::GetObject("NJettiness", key));
    REQUIRE(selector != nullptr);
    CHECK(selector->Trigger(particles) == (algorithm == "antikt"));
  }
}

TEST_CASE("NJettiness uses the Durham algorithm for kt in lepton collisions",
          "[PHASIC::NJettiness_Selector]")
{
  PHASIC_Test::SeedGlobals();
  ATOOLS::Particle_Info electron, gluon;
  electron.m_kfc = kf_e;
  gluon.m_kfc = kf_gluon;
  gluon.m_strong = 8;
  ATOOLS::Flavour flavour(gluon);
  const ATOOLS::Flavour e(electron);
  PHASIC_Test::Jet_Group jet_group(flavour);
  PHASIC_Test::Bunches bunches(e, e);
  PHASIC_Test::Process process({e, e, flavour, flavour, flavour}, 2);

  // Durham merges the soft gluon with its neighbour b at an angle of 0.3;
  // the two hard directions are a and b+soft.
  const auto gluon_at = [](double energy, double phi) {
    return ATOOLS::Vec4D(energy, energy * std::cos(phi), energy * std::sin(phi), 0.);
  };
  const ATOOLS::Vec4D a = gluon_at(50., 0.), b = gluon_at(45., 2.9),
                      soft = gluon_at(5., 2.6);
  ATOOLS::Selector_List particles(process.Flavours(),
    {ATOOLS::Vec4D(50., 0., 0., 50.), ATOOLS::Vec4D(50., 0., 0., -50.),
     a, b, soft}, 2);
  // For jet2 normalization, a two-parton axis contributes M^2/E.
  const double tau = (b + soft).Abs2() / (b + soft).E();
  for (const double cut : {tau / 2., 2. * tau}) {
    CAPTURE(cut);
    PHASIC::Selector_Key key;
    key.p_proc = &process;
    key.m_settings = ATOOLS::Scoped_Settings(
      "NJettiness: {N: 2, TauN: " + std::to_string(cut) + ", Algorithm: kt}");
    std::unique_ptr<PHASIC::Selector_Base> selector(
      PHASIC::Selector_Getter::GetObject("NJettiness", key));
    REQUIRE(selector != nullptr);
    bool trigger = false;
    REQUIRE_NOTHROW(trigger = selector->Trigger(particles));
    CHECK(trigger == (cut < tau));
  }
}
