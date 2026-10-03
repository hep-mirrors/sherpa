#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Test_Utils.H"
#include "PHASIC++/Scales/Tag_Setter.H"

namespace {
class Scale : public PHASIC::Scale_Setter_Base {
public:
  explicit Scale(const PHASIC::Scale_Setter_Arguments& args) : Scale_Setter_Base(args) {}
private:
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override
  {
    m_scale = {100., 200., 300.};
    return m_scale[0];
  }
};
}

TEST_CASE("Scale tags distinguish the massive transverse product from indexed scales",
          "[PHASIC::Tag_Setter]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle higgs(kf_h0);
  higgs.SetMass(4.);
  PHASIC_Test::Process process({ATOOLS::Flavour(kf_none), ATOOLS::Flavour(kf_none),
                                ATOOLS::Flavour(kf_h0), ATOOLS::Flavour(kf_h0),
                                ATOOLS::Flavour(kf_none)}, 2);
  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = &process;
  Scale setter(args);
  const ATOOLS::Vec4D_Vector momenta = {ATOOLS::Vec4D(10., 0., 0., 10.),
                                       ATOOLS::Vec4D(10., 0., 0., -10.),
                                       ATOOLS::Vec4D(5., 3., 0., 0.),
                                       ATOOLS::Vec4D(4., 0., 0., 0.),
                                       ATOOLS::Vec4D(2., 0., 2., 0.)};
  setter.CalculateScale(momenta);
  PHASIC::Tag_Setter tags(&setter);
  ATOOLS::Algebra_Interpreter calculator;
  calculator.SetTagReplacer(&tags);
  tags.SetCalculator(&calculator);
  tags.SetTags(&calculator);
  SECTION("P_TM2 is the squared geometric mean over massive final legs only") {
    calculator.Interprete("P_TM2");
    CHECK_THAT(calculator.Calculate()->Get<double>(),
               Catch::Matchers::WithinRel(momenta[2].MPerp() * momenta[3].MPerp(), 1.e-12));
  }
  SECTION("indexed scales and momentum tags retain their mapping") {
    for (size_t i = 0; i < setter.Scales().size(); ++i) {
      calculator.Interprete("MU_" + std::to_string(i) + "2");
      CHECK_THAT(calculator.Calculate()->Get<double>(),
                 Catch::Matchers::WithinRel(setter.Scales()[i], 1.e-12));
    }
    calculator.Interprete("PPerp2(p[4])");
    CHECK_THAT(calculator.Calculate()->Get<double>(),
               Catch::Matchers::WithinRel(momenta[4].PPerp2(), 1.e-12));
  }
  SECTION("the retired beam-thrust tag fails instead of returning P_SUM") {
    CHECK_THROWS([&] {
      calculator.Interprete("TAUB");
      calculator.Calculate();
    }());
  }
}

TEST_CASE("Hadronic H_T tags sum p_T or m_T of coloured final-state legs",
          "[PHASIC::Tag_Setter]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle electron(kf_e), gluon(kf_gluon), charm(kf_c);
  electron.Info().m_icharge = -3;
  gluon.Info().m_strong = 8;
  charm.SetMass(1.5);
  charm.Info().m_strong = 3;
  // e- g -> e- c cbar; the tags only read transverse quantities, so momentum
  // balance is irrelevant here
  PHASIC_Test::Process process({ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_gluon),
                                ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_c),
                                ATOOLS::Flavour(kf_c, true)}, 2);
  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = &process;
  Scale setter(args);
  const ATOOLS::Vec4D_Vector momenta = {ATOOLS::Vec4D(10., 0., 0., 10.),
                                       ATOOLS::Vec4D(10., 0., 0., -10.),
                                       ATOOLS::Vec4D(10., -1., 0., std::sqrt(99.)),
                                       ATOOLS::Vec4D(std::sqrt(7.25), 2., 0., 1.),
                                       ATOOLS::Vec4D(3.5, -1., 0., -3.)};
  setter.CalculateScale(momenta);
  PHASIC::Tag_Setter tags(&setter);
  ATOOLS::Algebra_Interpreter calculator;
  calculator.SetTagReplacer(&tags);
  tags.SetCalculator(&calculator);
  tags.SetTags(&calculator);
  const double htm = momenta[3].MPerp() + momenta[4].MPerp();
  calculator.Interprete("hH_T2");
  CHECK_THAT(calculator.Calculate()->Get<double>(), Catch::Matchers::WithinRel(9., 1.e-12));
  calculator.Interprete("hH_TM2");
  CHECK_THAT(calculator.Calculate()->Get<double>(),
             Catch::Matchers::WithinRel(htm * htm, 1.e-12));
}

TEST_CASE("Dressed scale tags reject arguments they cannot apply",
          "[PHASIC::Tag_Setter]")
{
  PHASIC_Test::SeedGlobals();
  PHASIC_Test::ScopedParticle electron(kf_e);
  electron.Info().m_icharge = -3;
  // e+e- -> e+e- has no photon to dress, so the dressing arguments used to be
  // ignored instead of validated
  PHASIC_Test::Process process({ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_e, true),
                                ATOOLS::Flavour(kf_e), ATOOLS::Flavour(kf_e, true)}, 2);
  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = &process;
  Scale setter(args);
  const ATOOLS::Vec4D_Vector momenta = {ATOOLS::Vec4D(50., 0., 0., 50.),
                                       ATOOLS::Vec4D(50., 0., 0., -50.),
                                       ATOOLS::Vec4D(50., 30., 0., 40.),
                                       ATOOLS::Vec4D(50., -30., 0., -40.)};
  setter.CalculateScale(momenta);
  PHASIC::Tag_Setter tags(&setter);
  ATOOLS::Algebra_Interpreter calculator;
  calculator.SetTagReplacer(&tags);
  tags.SetCalculator(&calculator);
  tags.SetTags(&calculator);
  for (const std::string function : {"DH_Tp2", "DH_Tln2"}) {
    // kt is advertised but not implemented; a missing dR is malformed
    for (const std::string arguments : {"(kt,0.1)", "(Cone)"}) {
      INFO(function + arguments);
      CHECK_THROWS([&] {
        calculator.Interprete(function + arguments);
        calculator.Calculate();
      }());
    }
  }
}
