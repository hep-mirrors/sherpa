#include <catch2/catch_all.hpp>

#include "Tests/Unit/PHASIC_Coupling_Test_Utils.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include <cmath>
#include <memory>

namespace {
class Couplings {
  PHASIC_Test::MasslessQCDParticles m_particles;
  MODEL::Running_AlphaS* p_previous_as;
  MODEL::Running_AlphaQED* p_previous_aqed;
  std::unique_ptr<MODEL::Running_AlphaS> m_as;
  MODEL::Running_AlphaQED m_aqed{1. / 137.};
public:
  Couplings()
  {
    PHASIC_Test::SeedGlobals();
    p_previous_as = MODEL::as;
    p_previous_aqed = MODEL::aqed;
    m_as.reset(new MODEL::Running_AlphaS(.118, 91.1876 * 91.1876, 0, 1));
    MODEL::as = m_as.get();
    MODEL::aqed = &m_aqed;
  }
  ~Couplings()
  {
    MODEL::as = p_previous_as;
    MODEL::aqed = p_previous_aqed;
    m_as.reset();
  }
};

PHASIC::Scale_Setter_Arguments ScaleArguments(PHASIC::Process_Base* process)
{
  PHASIC::Scale_Setter_Arguments args(nullptr);
  args.p_proc = process;
  return args;
}

class Scale : public PHASIC::Scale_Setter_Base {
public:
  explicit Scale(PHASIC::Process_Base* process)
    : Scale_Setter_Base(ScaleArguments(process))
  { m_scale = {100., 200., 300.}; }
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override { return m_scale[0]; }
};

class Process : public PHASIC_Test::Process {
public:
  Process() : PHASIC_Test::Process(ATOOLS::Flavour_Vector(4, ATOOLS::Flavour(kf_none)), 2)
  { p_scale = new Scale(this); }
};
}

TEST_CASE("Variable K-factor indexed and alias scales follow the current event or subevent",
          "[PHASIC::Variable_KFactor_Setter]")
{
  Couplings couplings;
  Process process;
  std::string expression;
  SECTION("indexed scale tags are initialized before parsing")
    expression = "VAR{MU_12+10*MU_02}";
  SECTION("scale aliases agree with indexed tags on each subevent")
    // Differences test alias consistency; the remaining terms distinguish
    // the selected scale vector rather than checking equality alone.
    expression = "VAR{MU_R2-MU_12+MU_F2-MU_02+MU_R2+10*MU_F2}";
  std::unique_ptr<PHASIC::KFactor_Setter_Base> setter(
    PHASIC::KFactor_Setter_Base::KFactor_Getter_Function::GetObject(
      "VAR", PHASIC::KFactor_Setter_Arguments(expression, &process)));
  REQUIRE(setter != nullptr);
  const auto expected = [](const std::vector<double>& scales) {
    return scales[ATOOLS::stp::ren] + 10. * scales[ATOOLS::stp::fac];
  };
  CHECK_THAT(setter->KFactor(), Catch::Matchers::WithinRel(expected(process.ScaleSetter()->Scales()), 1.e-12));
  ATOOLS::Vec4D_Vector momenta(4);
  ATOOLS::NLO_subevt subevent;
  subevent.m_n = momenta.size();
  subevent.p_mom = momenta.data();
  for (const auto& scales : {std::vector<double>{7., 11., 13.},
                            std::vector<double>{17., 23., 29.}}) {
    subevent.m_mu2 = scales;
    CHECK_THAT(setter->KFactor(subevent), Catch::Matchers::WithinRel(expected(scales), 1.e-12));
    CHECK_THAT(setter->LastKFactor(), Catch::Matchers::WithinRel(expected(scales), 1.e-12));
  }
  CHECK_THAT(setter->KFactor(), Catch::Matchers::WithinRel(expected(process.ScaleSetter()->Scales()), 1.e-12));
  const std::vector<double> parent{100., 200., 300.};
  for (size_t i = 0; i < parent.size(); ++i)
    CHECK_THAT(process.ScaleSetter()->Scales()[i], Catch::Matchers::WithinRel(parent[i], 1.e-12));
}

TEST_CASE("Variable K-factors stay finite on subevents that fail the cuts",
          "[PHASIC::Variable_KFactor_Setter]")
{
  Couplings couplings;
  Process process;
  std::string expression;
  SECTION("ratio of the couplings at the scale aliases") {
    expression = "VAR{Alpha_QCD(MU_R2)/Alpha_QCD(MU_F2)}";
  }
  SECTION("documented HEFT reweighting with an indexed scale") {
    expression = "VAR{sqr(Alpha_QCD(sqr(125))/Alpha_QCD(MU_12))}";
  }
  std::unique_ptr<PHASIC::KFactor_Setter_Base> setter(
    PHASIC::KFactor_Setter_Base::KFactor_Getter_Function::GetObject(
      "VAR", PHASIC::KFactor_Setter_Arguments(expression, &process)));
  REQUIRE(setter != nullptr);
  // a dipole subevent; Scale_Setter_Base gives it scales of -1 when it fails
  // the cuts, and its weight is zero
  ATOOLS::Vec4D_Vector momenta(4);
  ATOOLS::NLO_subevt dipole(momenta.size(), nullptr, nullptr, momenta.data(), 2, 3, 0);
  dipole.m_mu2.assign(dipole.m_mu2.size(), -1.);
  CHECK(std::isfinite(setter->KFactor(dipole)));
  // once it passes the cuts, the K-factor is evaluated at its scales
  dipole.m_trig = true;
  dipole.m_mu2 = {400., 900., 1600.};
  const double ren = MODEL::as->AlphaS(900.), fac = MODEL::as->AlphaS(400.);
  const double expected = expression.find("MU_R2") != std::string::npos
    ? ren / fac : std::pow(MODEL::as->AlphaS(125. * 125.) / ren, 2);
  CHECK_THAT(setter->KFactor(dipole), Catch::Matchers::WithinRel(expected, 1.e-12));
}
