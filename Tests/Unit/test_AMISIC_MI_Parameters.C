#include <catch2/catch_all.hpp>

#include "AMISIC++/Tools/MI_Parameters.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <cmath>

using namespace AMISIC;
using namespace ATOOLS;

namespace {
  // A bare unit-test binary never runs Sherpa's main, so the message singleton
  // is null and the global Run_Parameter is absent. MI_Parameters' ctor reads
  // rpa->gen.Ecms(), so create both once and pin the cms energy.
  void EnsureRuntime(double ecms)
  {
    if (!msg) msg = new Message();
    if (!rpa) rpa = new Run_Parameter();
    rpa->gen.SetEcms(ecms);
  }
}

TEST_CASE("MI_Parameters reads AMISIC settings and computes pt scales",
          "[AMISIC::MI_Parameters]") {

  // built-in proton defaults exercised below:
  //   pt_0(ref)=2.05, pt_0(IR)=0.5, pt_min(ref)=1.10, pt_min(IR)=1.0,
  //   E(ref)=7000, eta=0.026
  const double Sref = 7000. * 7000.;
  const double Scms = 13000. * 13000.;

  SECTION("CalculatePT02/PTmin2 scale with s and are floored at the IR value") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{"AMISIC: {}"});
    MI_Parameters p;

    // at the reference scale the regulator equals pt_(ref)^2
    CHECK(p.CalculatePT02(Sref)   == Catch::Approx(2.05 * 2.05));
    CHECK(p.CalculatePTmin2(Sref) == Catch::Approx(1.10 * 1.10));
    // far below the reference scale the pt_min^2 IR floor (Max) takes over
    CHECK(p.CalculatePTmin2(1.0)  == Catch::Approx(1.0));
    // pt_0^2 still scales here: its IR floor 0.5^2 lies far below the scaled
    // value for the default eta, so the Max does not clamp it
    CHECK(p.CalculatePT02(1.0)
          == Catch::Approx(2.05 * 2.05 * std::pow(1.0 / Sref, 2 * 0.026)));
    // the scaling exponent is 2*eta
    CHECK(p.CalculatePT02(4. * Sref) / p.CalculatePT02(Sref)
          == Catch::Approx(std::pow(4.0, 2 * 0.026)));
    // s<0 falls back to the run cms energy squared
    CHECK(p.CalculatePT02(-1.) == Catch::Approx(p.CalculatePT02(Scms)));
  }

  SECTION("pt_min/pt_0 defaults are the energy-scaled computed values") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{"AMISIC: {}"});
    MI_Parameters p;
    CHECK(p("pt_min") == Catch::Approx(std::sqrt(p.CalculatePTmin2(Scms))));
    CHECK(p("pt_0")   == Catch::Approx(std::sqrt(p.CalculatePT02(Scms))));
  }

  SECTION("unset parameters take their built-in defaults") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{"AMISIC: {}"});
    MI_Parameters p;
    CHECK(p("pt_0(ref)")    == Catch::Approx(2.05));
    CHECK(p("eta")          == Catch::Approx(0.026));
    CHECK(p("SigmaND_Norm") == Catch::Approx(0.44));
    CHECK(p["nMC_points"]   == 100000);
    CHECK(p["nS_bins"]      == 40);
  }

  SECTION("enum and flag settings are parsed") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{
        "AMISIC: {EVENT_TYPE: Elastic, MU_R_SCHEME: PT_with_Raps, "
        "MU_F_SCHEME: PT, TwoPionInterference: 3}"});
    MI_Parameters p;
    CHECK(p.GetEvtType()          == evt_type::Elastic);
    CHECK(p.GetScaleRScheme()     == scale_scheme::PT_with_Raps);
    CHECK(p.GetScaleFScheme()     == scale_scheme::PT);
    CHECK(p.GetTwoPionTreatment() == two_pions::rho_omega_cont);
  }

  SECTION("EVENT_TYPE accepts MinimumBias and Non-Perturbative") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{"AMISIC: {EVENT_TYPE: MinimumBias}"});
    CHECK(MI_Parameters{}.GetEvtType() == evt_type::AllMinimumBias);
    Settings::InitializeMainSettings(std::string{"AMISIC: {EVENT_TYPE: Non-Perturbative}"});
    CHECK(MI_Parameters{}.GetEvtType() == evt_type::NonPerturbative);
  }

  SECTION("invalid input is rejected") {
    EnsureRuntime(13000.);
    Settings::InitializeMainSettings(std::string{"AMISIC: {TRIGGER: [1, 2, 3]}"});
    REQUIRE_THROWS(MI_Parameters{});               // > 2 trigger flavours

    Settings::InitializeMainSettings(std::string{"AMISIC: {}"});
    MI_Parameters p;
    REQUIRE_THROWS(p("no_such_parameter"));
    REQUIRE_THROWS(p["no_such_flag"]);
  }
}
