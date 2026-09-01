#include <catch2/catch_all.hpp>

#include <sstream>

#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace REMNANTS;
using namespace ATOOLS;

namespace {
  // The bare unit-test binary never runs Sherpa's main, so the global message
  // singleton is null; the settings machinery dereferences it. Create it once.
  void EnsureMessage()
  {
    if (!msg) msg = new Message();
  }

  // Register the handful of particles that Remnants_Parameters builds in its
  // constructor (and a few exotics used by the fallback specs). The KF table is
  // a process-global shared with every other test case in this binary, and
  // AddParticle throws on a duplicate kf code, so skip codes already present.
  // AddParticle reads PARTICLE_DATA from the main settings, so the settings
  // singleton must be seeded first.
  void Add(kf_code kfc, double mass, int icharge, int spin,
           const std::string& name, const std::string& antiname)
  {
    if (s_kftable.find(kfc) != s_kftable.end()) return;
    AddParticle(kfc, mass, 0., 0., icharge, 0, spin, 0, 1, 1, 1,
                name, antiname, name, antiname, 0, 0);
  }

  void EnsureParticles()
  {
    //  kfc   mass    ich spin  names
    Add(2212, 0.938,   3,  1, "P+", "P-");
    Add(2112, 0.940,   0,  1, "n", "nb");
    Add( 211, 0.140,   3,  0, "pi+", "pi-");
    Add( 111, 0.135,   0,  0, "pi", "pi");
    Add(  22, 0.,      0,  2, "P", "P");
    Add(  11, 5.1e-4, -3,  1, "e-", "e+");
    Add(  13, 0.106,  -3,  1, "mu-", "mu+");
    Add(3122, 1.116,   0,  1, "Lambda", "Lambdab");
    Add( 321, 0.494,   3,  0, "K+", "K-");
  }
}

TEST_CASE("Remnants_Parameters serves the hardcoded per-flavour default tune",
          "[REMNANTS::Remnants_Parameters]") {

  SECTION("no YAML override -> the built-in defaults per flavour") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{"REMNANTS: {}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus), pion(kf_pi_plus), photon(kf_photon),
                  electron(kf_e);

    // proton (nucleon tune)
    CHECK(rp.KT_Form(proton)    == primkT_form::gauss_limited);
    CHECK(rp.KT_Recoil(proton)  == primkT_recoil::beam_vs_shower);
    CHECK(rp.Matter_Form(proton) == matter_form::double_gaussian);
    CHECK(rp.Get(proton, "MATTER_RADIUS_1") == Catch::Approx(0.85));
    CHECK(rp.Get(proton, "SOFT_EXPONENT")   == Catch::Approx(0.08));

    // charged pion (meson tune)
    CHECK(rp.Matter_Form(pion)  == matter_form::single_gaussian);
    CHECK(rp.Get(pion, "MATTER_RADIUS_1") == Catch::Approx(0.75));
    CHECK(rp.Get(pion, "SOFT_EXPONENT")   == Catch::Approx(0.0).margin(1e-12));

    // photon tune
    CHECK(rp.Matter_Form(photon) == matter_form::single_gaussian);
    CHECK(rp.Get(photon, "MATTER_RADIUS_1") == Catch::Approx(0.75));

    // electron (lepton tune)
    CHECK(rp.KT_Form(electron)    == primkT_form::none);
    CHECK(rp.Matter_Form(electron) == matter_form::none);
    CHECK(rp.Get(electron, "MATTER_RADIUS_1") == Catch::Approx(1e-12));
  }

  SECTION("flavours absent from the table fall back by hadron class") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{"REMNANTS: {}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour lambda(3122), kaon(321), none(kf_none);

    // an exotic baryon (Lambda) is not in the tune -> proton defaults
    CHECK(rp.KT_Form(lambda)     == primkT_form::gauss_limited);
    CHECK(rp.Matter_Form(lambda) == matter_form::double_gaussian);
    CHECK(rp.Get(lambda, "MATTER_RADIUS_1") == Catch::Approx(0.85));

    // an exotic meson (kaon) is not in the tune -> pion defaults
    CHECK(rp.Matter_Form(kaon) == matter_form::single_gaussian);
    CHECK(rp.Get(kaon, "MATTER_RADIUS_1") == Catch::Approx(0.75));

    // the empty flavour: enum getters short-circuit to none (production never
    // calls the scalar Get() with kf_none, which has no registered kf here)
    CHECK(rp.KT_Form(none)     == primkT_form::none);
    CHECK(rp.KT_Recoil(none)   == primkT_recoil::none);
    CHECK(rp.Matter_Form(none) == matter_form::none);
  }

  SECTION("a per-pid YAML entry overrides the default for that flavour") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {MATTER_RADIUS_1: 1.234, MATTER_FORM: Single_Gaussian}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus);

    // overridden keys take the YAML value ...
    CHECK(rp.Get(proton, "MATTER_RADIUS_1") == Catch::Approx(1.234));
    CHECK(rp.Matter_Form(proton) == matter_form::single_gaussian);
    // ... while un-overridden keys keep their defaults
    CHECK(rp.Get(proton, "SOFT_EXPONENT") == Catch::Approx(0.08));
    CHECK(rp.KT_Form(proton) == primkT_form::gauss_limited);
  }

  SECTION("a per-pid override now applies to the antiparticle too") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {MATTER_RADIUS_1: 1.234}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour antiproton = Flavour(kf_p_plus).Bar();

    // keyed on |kfcode|, so the 2212 override reaches the antiproton
    CHECK(rp.Get(antiproton, "MATTER_RADIUS_1") == Catch::Approx(1.234));
  }

  SECTION("SOFT_EXPONENT is read from YAML") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {SOFT_EXPONENT: 0.5}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    CHECK(rp.Get(Flavour(kf_p_plus), "SOFT_EXPONENT") == Catch::Approx(0.5));
  }
}

TEST_CASE("Remnants_Parameters serves the on-the-fly variation vectors",
          "[REMNANTS::Remnants_Parameters]") {

  // The four matter-distribution parameters may be given as a YAML list, with
  // the nominal value first and the remaining entries used as on-the-fly
  // reweighting variations (consumed by Form_Factor). Get() must keep
  // returning the nominal value for all of them.

  SECTION("no YAML override -> the single-element default vector") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{"REMNANTS: {}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus);

    for (const auto& key : {std::string("MATTER_RADIUS_1"),
                            std::string("MATTER_RADIUS_2"),
                            std::string("MATTER_FRACTION_1"),
                            std::string("SOFT_EXPONENT")}) {
      const std::vector<double> var = rp.GetVariationVector(proton, key);
      REQUIRE(var.size() == 1);
      // the nominal entry is exactly what the scalar getter reports
      CHECK(var.front() == Catch::Approx(rp.Get(proton, key)));
    }

    CHECK(rp.GetVariationVector(proton, "MATTER_RADIUS_1").front()
          == Catch::Approx(0.85));
    CHECK(rp.GetVariationVector(proton, "MATTER_FRACTION_1").front()
          == Catch::Approx(0.65));
  }

  SECTION("a YAML list yields the full variation vector, Get() the nominal") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {MATTER_RADIUS_1: [0.85, 0.80, 0.90], "
        "MATTER_FRACTION_1: [0.65, 0.60, 0.70]}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus);

    const std::vector<double> radius = rp.GetVariationVector(proton, "MATTER_RADIUS_1");
    REQUIRE(radius.size() == 3);
    CHECK(radius[0] == Catch::Approx(0.85));
    CHECK(radius[1] == Catch::Approx(0.80));
    CHECK(radius[2] == Catch::Approx(0.90));

    const std::vector<double> fraction = rp.GetVariationVector(proton, "MATTER_FRACTION_1");
    REQUIRE(fraction.size() == 3);
    CHECK(fraction[0] == Catch::Approx(0.65));

    // the scalar getter returns the nominal (first) entry, not the last read
    CHECK(rp.Get(proton, "MATTER_RADIUS_1")   == Catch::Approx(0.85));
    CHECK(rp.Get(proton, "MATTER_FRACTION_1") == Catch::Approx(0.65));

    // an un-listed variation key keeps its single-element default
    CHECK(rp.GetVariationVector(proton, "SOFT_EXPONENT").size() == 1);
  }

  SECTION("a scalar YAML value is auto-wrapped into a one-element vector") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {MATTER_RADIUS_1: 1.234}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus);

    const std::vector<double> radius = rp.GetVariationVector(proton, "MATTER_RADIUS_1");
    REQUIRE(radius.size() == 1);
    CHECK(radius.front() == Catch::Approx(1.234));
    CHECK(rp.Get(proton, "MATTER_RADIUS_1") == Catch::Approx(1.234));
  }

  SECTION("keys without variations report an empty vector") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{"REMNANTS: {}"});
    EnsureParticles();
    Remnants_Parameters rp;

    const Flavour proton(kf_p_plus);

    // only the four matter-distribution parameters carry variations; the kT
    // parameters are plain scalars and must not be reported as varied.
    CHECK(rp.GetVariationVector(proton, "SHOWER_INITIATOR_MEAN").empty());
    CHECK(rp.GetVariationVector(proton, "BEAM_SPECTATOR_SIGMA").empty());
    CHECK(rp.GetVariationVector(proton, "REFERENCE_ENERGY").empty());
    // ... and neither does a keyword that is not a parameter at all
    CHECK(rp.GetVariationVector(proton, "NO_SUCH_KEY").empty());
  }

  SECTION("variation lists follow the hadron-class fallback and |kfcode|") {
    EnsureMessage();
    Settings::InitializeMainSettings(std::string{
        "REMNANTS: {2212: {MATTER_RADIUS_1: [0.85, 0.80]}}"});
    EnsureParticles();
    Remnants_Parameters rp;

    // keyed on |kfcode|, so the 2212 list reaches the antiproton
    const std::vector<double> anti =
        rp.GetVariationVector(Flavour(kf_p_plus).Bar(), "MATTER_RADIUS_1");
    REQUIRE(anti.size() == 2);
    CHECK(anti[1] == Catch::Approx(0.80));

    // an exotic baryon is not in the tune and has no 3122 YAML entry, so it
    // falls back to the proton *default*, not the proton's YAML override
    const std::vector<double> lambda =
        rp.GetVariationVector(Flavour(3122), "MATTER_RADIUS_1");
    REQUIRE(lambda.size() == 1);
    CHECK(lambda.front() == Catch::Approx(0.85));
  }
}

namespace {
  // Stream an enum out and check operator<< emits a single whitespace-free
  // token (ToString truncates at the first whitespace), then parse it back
  // with operator>> and require the exact original value. This is the
  // round-trip contract every settings enum must satisfy.
  template <typename Enum>
  Enum RoundTrip(Enum value)
  {
    std::stringstream os;
    os << value;
    std::string tok;
    std::stringstream(os.str()) >> tok;
    CHECK(tok == os.str());        // single token, no embedded whitespace
    Enum back;
    std::stringstream(os.str()) >> back;
    return back;
  }
}

TEST_CASE("REMNANTS settings enums round-trip through their stream operators",
          "[REMNANTS::Remnants_Parameters]") {
  // operator<< / operator>> must be exact inverses for every enumerator that
  // can appear as a default or a YAML value, so that an enum default written
  // via ToString and parsed back through operator>> never silently collapses
  // to a wrong value (this regression guards the round-trip contract that the
  // matter_form parser previously violated).

  SECTION("primkT_form") {
    for (auto v : {primkT_form::none, primkT_form::gauss,
                   primkT_form::gauss_limited, primkT_form::dipole,
                   primkT_form::dipole_limited})
      CHECK(RoundTrip(v) == v);
  }

  SECTION("primkT_recoil") {
    // primkT_recoil::none is never written as a default (returned directly for
    // kf_none) and intentionally has no operator<< token, so it is excluded.
    for (auto v : {primkT_recoil::democratic, primkT_recoil::beam_vs_shower})
      CHECK(RoundTrip(v) == v);
  }

  SECTION("matter_form") {
    for (auto v : {matter_form::none, matter_form::single_gaussian,
                   matter_form::double_gaussian,
                   matter_form::x_dependent_gaussian, matter_form::unknown})
      CHECK(RoundTrip(v) == v);
  }
}
