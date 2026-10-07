#include <catch2/catch_all.hpp>

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include <cmath>
#include <memory>
#include <sstream>
#include <string>

using namespace ATOOLS;
using namespace MODEL;

namespace {

  const double alpha0 = 1. / 137.03599976;
  const double mz2    = 91.1876 * 91.1876;

  // The unit tests share one process, and other tests register or clear
  // particles, so make sure the leptons and the top exist with the masses
  // used here every time instead of only once.
  void EnsureParticle(kf_code kf, double mass, int icharge, int strong,
                      const char *name, const char *anti)
  {
    if (s_kftable.find(kf) == s_kftable.end())
      AddParticle(kf, mass, .0, .0, icharge, strong, 1, 0, 1, 1, 0, name, anti,
                  name, anti);
    Flavour(kf).SetMass(mass);
  }

  void Boot()
  {
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Settings::InitializeMainSettings("");
    EnsureParticle(kf_e, 0.000511, -3, 0, "e-", "e+");
    EnsureParticle(kf_mu, 0.105, -3, 0, "mu-", "mu+");
    EnsureParticle(kf_tau, 1.777, -3, 0, "tau-", "tau+");
    EnsureParticle(kf_t, 173.2, 2, 1, "t", "tb");
  }

  // alpha_QED(t) for a given Alpha_QED: VPMODE setting
  double AlphaQED(const std::string &mode, double t)
  {
    Boot();
    Settings::InitializeMainSettings("Alpha_QED: {VPMODE: " + mode + "}");
    EnsureParticle(kf_t, 173.2, 2, 1, "t", "tb");
    Running_AlphaQED aqed(alpha0);
    return aqed(t);
  }

  vpmode::code Parse(const std::string &tag)
  {
    std::istringstream str(tag);
    vpmode::code mode;
    str >> mode;
    return mode;
  }

  double Delta(const std::string &mode, double t)
  {
    return 1. - alpha0 / AlphaQED(mode, t);
  }

}  // namespace

TEST_CASE("VPMODE tags are parsed, and unknown ones are rejected",
          "[MODEL::Running_AlphaQED]")
{
  Boot();
  CHECK(Parse("Default") == vpmode::builtin);
  CHECK(Parse("None") == vpmode::off);
  CHECK(Parse("Full") == vpmode::full);
  CHECK(Parse("HP") == vpmode::hp);
  CHECK(Parse("LP") == vpmode::lp);
  // tags are case sensitive, like all other Sherpa settings
  CHECK_THROWS(Parse("full"));
  CHECK_THROWS(Parse("hp"));
  CHECK_THROWS(Parse("Fuul"));
  CHECK_THROWS(Parse("Legacy"));
}

TEST_CASE("Default VPMODE keeps the built-in running of alpha_QED",
          "[MODEL::Running_AlphaQED]")
{
  Boot();
  Settings::InitializeMainSettings("");
  Running_AlphaQED aqed(alpha0);
  CHECK(aqed.m_mode == vpmode::builtin);

  SECTION("alpha(0) is alpha0")
  {
    CHECK(aqed(0.) == alpha0);
  }

  SECTION("the value at m_Z^2 is the documented 1/128.802")
  {
    CHECK(1. / aqed(mz2) == Catch::Approx(128.802).epsilon(1.e-4));
  }

  SECTION("it depends on |t| only")
  {
    CHECK(aqed(-mz2) == Catch::Approx(aqed(mz2)).epsilon(1.e-14));
  }
}

TEST_CASE("VPMODE None freezes alpha_QED at alpha(0)",
          "[MODEL::Running_AlphaQED]")
{
  for (double t : {-1.e4, -1., 0., 1., 8315.18})
    CHECK(AlphaQED("None", t) == alpha0);
}

#ifdef USING__HADALPHAQED

TEST_CASE("Jegerlehner hadr5x modes at the Z pole",
          "[MODEL::Running_AlphaQED]")
{
  const double dhad = Delta("HP", mz2), dlep = Delta("LP", mz2),
               dfull = Delta("Full", mz2);

  SECTION("hadronic part is Delta alpha_had^(5)(m_Z) = 0.0277")
  {
    // 0.02766 +- 0.00007 (Jegerlehner); passing Q^2 instead of the energy to
    // hadr5x gave 0.0542, so the tolerance is tight enough to catch that
    CHECK(dhad == Catch::Approx(0.02766).margin(2.e-4));
  }

  SECTION("leptonic part is Delta alpha_lep(m_Z) = 0.0315")
  {
    CHECK(dlep == Catch::Approx(0.03150).margin(2.e-4));
  }

  SECTION("Full is leptons + hadrons + top")
  {
    // the top contribution at m_Z is about -7e-5
    CHECK(dfull == Catch::Approx(dlep + dhad).margin(2.e-4));
    CHECK(dfull > dhad);
  }

  SECTION("Full and the built-in parameterisation agree to a few 1e-4")
  {
    CHECK(AlphaQED("Full", mz2) ==
          Catch::Approx(AlphaQED("Default", mz2)).epsilon(5.e-3));
  }

  SECTION("spacelike and timelike arguments differ (hadr5x is signed)")
  {
    CHECK(AlphaQED("HP", mz2) != AlphaQED("HP", -mz2));
  }
}

#else

TEST_CASE("Jegerlehner modes need SHERPA_ENABLE_ALPHAQEDHAD",
          "[MODEL::Running_AlphaQED]")
{
  Boot();
  Settings::InitializeMainSettings("Alpha_QED: {VPMODE: Full}");
  CHECK_THROWS(Running_AlphaQED(alpha0));
}

#endif
