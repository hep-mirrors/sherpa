#include <catch2/catch_all.hpp>
#include "Tests/Unit/PHASIC_Selector_Test_Utils.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Math/Function_Base.H"
#include <memory>

namespace {
ATOOLS::Vec4D_Vector recorded_momenta;

class AlphaQED : public ATOOLS::Function_Base {
public:
  AlphaQED() { m_name = "Alpha_QED"; m_defval = 1. / 137.; }
};

class Model : public MODEL::Model_Base {
  MODEL::Model_Base* p_previous;
public:
  Model() : Model_Base(false), p_previous(MODEL::s_model)
  {
    Constants()->emplace("alpha_QED", 1. / 137.);
    p_functions->emplace("alpha_QED", new AlphaQED());
    MODEL::s_model = this;
  }
  ~Model() override { MODEL::s_model = p_previous; }
private:
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

class RecordingLoop : public PHASIC::Virtual_ME2_Base {
public:
  explicit RecordingLoop(const PHASIC::Process_Info& info)
    : Virtual_ME2_Base(info, info.ExtractFlavours()) {}
  void Calc(const ATOOLS::Vec4D_Vector& momenta) override
  {
    recorded_momenta = momenta;
    m_born = 2.;
    m_res.Finite() = 3.;
  }
};

class Scale : public PHASIC::Scale_Setter_Base {
public:
  explicit Scale(const PHASIC::Scale_Setter_Arguments& args) : Scale_Setter_Base(args)
  {
    m_scale.assign(ATOOLS::stp::size, 100.);
    m_fac.assign(2, 1.);
  }
private:
  double Calculate(const ATOOLS::Vec4D_Vector&, const size_t&) override { return 100.; }
};

class Process : public PHASIC_Test::Process {
public:
  Process() : PHASIC_Test::Process(ATOOLS::Flavour_Vector(4, ATOOLS::Flavour(kf_none)), 2)
  {
    Info().m_loopgenerator = "PHASIC_EWVirt_RecordingLoop";
    for (size_t i = 0; i < 2; ++i) {
      Info().m_ii.m_ps.emplace_back(ATOOLS::Flavour(kf_none));
      Info().m_fi.m_ps.emplace_back(ATOOLS::Flavour(kf_none));
    }
    SetMaxOrders({0., 0.});
    SetMinOrders({0., 0.});
    PHASIC::Scale_Setter_Arguments args(MODEL::s_model);
    args.p_proc = this;
    p_scale = new Scale(args);
  }
};
}

DECLARE_VIRTUALME2_GETTER(RecordingLoop, "!PHASIC_EWVirt_RecordingLoop");

PHASIC::Virtual_ME2_Base* ATOOLS::Getter<PHASIC::Virtual_ME2_Base,
  PHASIC::Process_Info, RecordingLoop>::operator()(const PHASIC::Process_Info& info) const
{
  return info.m_loopgenerator == "PHASIC_EWVirt_RecordingLoop" ? new RecordingLoop(info) : nullptr;
}

TEST_CASE("EWVirt forwards current subevent momenta to its loop provider",
          "[PHASIC::EWVirtKFactor_Setter]")
{
  PHASIC_Test::SeedGlobals();
  Model model;
  Process process;
  std::unique_ptr<PHASIC::KFactor_Setter_Base> setter(
    PHASIC::KFactor_Setter_Base::KFactor_Getter_Function::GetObject(
      "EWVirt", PHASIC::KFactor_Setter_Arguments("EWVirt", &process)));
  REQUIRE(setter != nullptr);
  ATOOLS::Vec4D_Vector momenta = {
    {-10., 0., 0., -10.}, {-10., 0., 0., 10.},
    {10., 6., 0., 8.}, {10., -6., 0., -8.}};
  ATOOLS::NLO_subevt subevent;
  subevent.m_n = momenta.size();
  subevent.p_mom = momenta.data();
  for (size_t event = 0; event < 2; ++event) {
    CAPTURE(event);
    CHECK_THAT(setter->KFactor(subevent),
               Catch::Matchers::WithinRel(1. + 3. / (137. * 2. * M_PI), 1.e-12));
    REQUIRE(recorded_momenta.size() == momenta.size());
    for (size_t i = 0; i < momenta.size(); ++i)
      CHECK(recorded_momenta[i] == (i < process.NIn() ? -momenta[i] : momenta[i]));
    CHECK(subevent.p_mom[0] == momenta[0]);
    for (auto& momentum : momenta) momentum *= 2.;
  }
  recorded_momenta.clear();
}
