#include "BEAM/Spectra/EPA.H"

#include <catch2/catch_all.hpp>

#include "Test_Setup.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "BEAM/Main/Collider_Weight.H"
#include "BEAM/Main/Transverse_Profile.H"
#include "BEAM/Spectra/Monochromatic.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <vector>

using namespace ATOOLS;
using namespace BEAM;

namespace {

  using Type = Transverse_Profile_Type;

  // Configuration cases deliberately initialise their own settings. Clear
  // their YAML on exit so no profile overrides leak into other test files.
  struct BeamSettings {
    explicit BeamSettings(const std::string& yaml)
    {
      Test_Setup::InitGlobals();
      Settings::InitializeMainSettings(yaml);
      Test_Setup::RegisterParticles();
    }
    ~BeamSettings() { Settings::InitializeMainSettings(""); }
  };

  // A mathematical test emitter: its normalised Gaussian is NOT an EPA flux.
  // Marking the spectrum EPA exercises the real overlap path and verifies
  // that Beam_Base does not replace a specialised emission sampler.
  class GaussianEmitter : public Beam_Base {
    Transverse_Profile m_emission;
  public:
    GaussianEmitter(Flavour beam, double rms_fm, int dir)
        : Beam_Base(beamspectrum::EPA, beam, 100., 0., dir),
          m_emission(Type::Gaussian, rms_fm / rpa->hBar_c()) {}
    Beam_Base* Copy() override { return new GaussianEmitter(*this); }
    bool CalculateWeight(double, double) override { return true; }
    double ImpactParameter() const override
    {
      return m_emission.SampleRadius(m_eran);
    }
  };

  class TestKinematics : public Kinematics_Base {
    void InitIntegration() override {}
  public:
    TestKinematics(Beam_Base* a, Beam_Base* b) : Kinematics_Base({a, b}) {}
    bool operator()(Vec4D_Vector&) override { return true; }
    void AssignKeys(Integration_Info*) override {}
    void SetLimits() override {}
  };

  // Deterministic low-discrepancy integration, independent of Sherpa's RNG.
  double RadicalInverse(unsigned int index, unsigned int base)
  {
    double result = 0., place = 1.;
    while (index) {
      place /= base;
      result += (index % base) * place;
      index /= base;
    }
    return result;
  }

} // namespace

TEST_CASE("Transverse profile models round-trip through settings",
          "[BEAM::Transverse_Profile]")
{
  Test_Setup::InitGlobals();
  for (Type type : {Type::Gaussian, Type::UniformSphere, Type::WoodsSaxon}) {
    std::ostringstream os;
    os << type;
    CHECK(os.str().find_first_of(" \t\n") == std::string::npos);
    Type parsed = Type::Gaussian;
    std::istringstream input(os.str());
    input >> parsed;
    CHECK_FALSE(input.fail());
    CHECK(parsed == type);

    // Use a local Settings instance: this check does not replace the singleton.
    Settings settings("");
    CHECK(settings["TYPE"].SetDefault(type).Get<Type>() == type);
  }
  Type type = Type::Gaussian;
  std::istringstream bad("UnknownProfile");
  bad >> type;
  CHECK(bad.fail());
}

TEST_CASE("Transverse profile inverse CDFs have the intended support and quantiles",
          "[BEAM::Transverse_Profile]")
{
  const double rms = 2.3;
  for (Type type : {Type::Gaussian, Type::UniformSphere}) {
    Transverse_Profile profile(type, rms);
    CHECK(profile.Type() == type);
    CHECK(profile.RMSRadius() == rms);
    CHECK(profile.SampleRadius(0.) == 0.);
    double previous = 0.;
    for (double u : {1.e-6, .01, .25, .5, .9, .999999}) {
      const double s = profile.SampleRadius(u);
      CHECK(s > previous);
      previous = s;
      const double cdf = type == Type::Gaussian
          ? -std::expm1(-s * s / (rms * rms))
          : 1. - std::pow(1. - s * s / (2.5 * rms * rms), 1.5);
      CHECK(cdf == Catch::Approx(u).epsilon(1.e-9));
    }
    const double edge = profile.SampleRadius(std::nextafter(1., 0.));
    CHECK(std::isfinite(edge));
    if (type == Type::UniformSphere)
      CHECK(edge <= std::sqrt(2.5) * rms);
    else
      CHECK(edge > 6. * rms); // no arbitrary Gaussian tail truncation

    const double tiny = 1.e-24;
    const double expected = rms * std::sqrt(tiny *
        (type == Type::Gaussian ? 1. : 5. / 3.));
    CHECK(profile.SampleRadius(tiny) == Catch::Approx(expected).epsilon(1.e-12));
    CHECK(profile.SampleRadius(.37) == profile.SampleRadius(.37));
  }
}

TEST_CASE("Transverse profile sampling includes the radial area measure",
          "[BEAM::Transverse_Profile]")
{
  const double rms = 1.7;
  const unsigned int n = 200000;
  for (Type type : {Type::Gaussian, Type::UniformSphere}) {
    Transverse_Profile profile(type, rms);
    double second = 0., fourth = 0.;
    for (unsigned int i = 0; i < n; ++i) {
      const double s = profile.SampleRadius((i + .5) / n);
      second += s * s;
      fourth += s * s * s * s;
    }
    CHECK(second / n == Catch::Approx(rms * rms).epsilon(1.e-5));
    const double expected = (type == Type::Gaussian ? 2. : 10. / 7.) *
                            std::pow(rms, 4);
    CHECK(fourth / n == Catch::Approx(expected).epsilon(1.e-4));
  }
}

TEST_CASE("Transverse profiles handle pointlike beams and reject invalid inputs",
          "[BEAM::Transverse_Profile]")
{
  Test_Setup::InitGlobals();
  const double inf = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (Type type : {Type::Gaussian, Type::UniformSphere}) {
    Transverse_Profile point(type, 0.);
    for (double u : {0., .5, std::nextafter(1., 0.)})
      CHECK(point.SampleRadius(u) == 0.);
    for (double bad : {-1., inf, nan})
      CHECK_THROWS(Transverse_Profile(type, bad));
    Transverse_Profile profile(type, 1.);
    for (double bad : {-1.e-9, 1., 1.1, inf, nan})
      CHECK_THROWS(profile.SampleRadius(bad));
  }
  CHECK_THROWS(Transverse_Profile(static_cast<Type>(99), 1.));
}

TEST_CASE("Beam profiles derive RMS sizes from beam radii without changing weights",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings("");
  Monochromatic proton(Flavour(kf_p_plus), 100., 0., 1);
  Monochromatic electron(Flavour(kf_e), 100., 0., -1);
  const double expected = std::sqrt(2. / 3.) * proton.Beam().Radius() / rpa->hBar_c();
  CHECK(proton.TransverseProfile().Type() == Type::Gaussian);
  CHECK(proton.TransverseProfile().RMSRadius() == Catch::Approx(expected));
  proton.SetERan(.5);
  CHECK(proton.ImpactParameter() == Catch::Approx(expected * std::sqrt(std::log(2.))));
  electron.SetERan(.7);
  CHECK(electron.ImpactParameter() == 0.);
  CHECK(proton.CalculateWeight(1., 1.));
  CHECK(proton.Weight() == 1.);
  CHECK(electron.Weight() == 1.);
}

TEST_CASE("Beam profile settings are per beam and survive spectrum copies",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_TRANSVERSE_PROFILE:\n"
      "  TYPE: [UniformSphere, Gaussian]\n"
      "  RMS_RADIUS: [2.0, 3.0]\n");
  const double original_radius = Flavour(kf_p_plus).Radius();
  Monochromatic first(Flavour(kf_p_plus), 100., 0., 1);
  Monochromatic second(Flavour(kf_p_plus), 100., 0., -1);
  CHECK(first.TransverseProfile().Type() == Type::UniformSphere);
  CHECK(second.TransverseProfile().Type() == Type::Gaussian);
  CHECK(first.TransverseProfile().RMSRadius() == Catch::Approx(2. / rpa->hBar_c()));
  CHECK(second.TransverseProfile().RMSRadius() == Catch::Approx(3. / rpa->hBar_c()));
  CHECK(Flavour(kf_p_plus).Radius() == original_radius);
  for (Beam_Base* beam : {static_cast<Beam_Base*>(&first),
                          static_cast<Beam_Base*>(&second)}) {
    std::unique_ptr<Beam_Base> copy(beam->Copy());
    CHECK(copy->TransverseProfile().Type() == beam->TransverseProfile().Type());
    CHECK(copy->TransverseProfile().RMSRadius() == beam->TransverseProfile().RMSRadius());
    beam->SetERan(.4);
    copy->SetERan(.4);
    CHECK(copy->ImpactParameter() == beam->ImpactParameter());
  }
}

TEST_CASE("Beam profile scalar settings broadcast and explicit zero is pointlike",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_TRANSVERSE_PROFILE: {TYPE: UniformSphere, RMS_RADIUS: 0}\n");
  for (int dir : {1, -1}) {
    Monochromatic beam(Flavour(kf_p_plus), 100., 0., dir);
    CHECK(beam.TransverseProfile().Type() == Type::UniformSphere);
    CHECK(beam.TransverseProfile().RMSRadius() == 0.);
    beam.SetERan(.9);
    CHECK(beam.ImpactParameter() == 0.);
  }
}

TEST_CASE("Programmatic beam profile overrides are honoured without a YAML entry",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings("");
  Settings::GetMainSettings()["BEAM_TRANSVERSE_PROFILE"]["RMS_RADIUS"]
      .OverrideVector(std::vector<double>{.6, 1.2});
  for (int dir : {1, -1}) {
    Monochromatic beam(Flavour(kf_p_plus), 100., 0., dir);
    CHECK(beam.TransverseProfile().RMSRadius() ==
          Catch::Approx((dir > 0 ? .6 : 1.2) / rpa->hBar_c()));
  }
}

TEST_CASE("Invalid beam profile configuration fails during setup",
          "[BEAM::Transverse_Profile]")
{
  for (const std::string& yaml : {
      "BEAM_TRANSVERSE_PROFILE: {TYPE: UnknownProfile}",
      "BEAM_TRANSVERSE_PROFILE: {TYPE: [Gaussian, Gaussian, Gaussian]}",
      "BEAM_TRANSVERSE_PROFILE: {RMS_RADIUS: -1}",
      "BEAM_TRANSVERSE_PROFILE: {RMS_RADIUS: nan}",
      "BEAM_TRANSVERSE_PROFILE: {RMS_RADIUS: []}",
      "BEAM_TRANSVERSE_PROFILE: {RMS_RADIUS: [1, 2, 3]}"}) {
    BeamSettings settings(yaml);
    CHECK_THROWS(Monochromatic(Flavour(kf_p_plus), 100., 0., 1));
  }
}

TEST_CASE("Specialised emission profiles are not replaced by generic beam settings",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_TRANSVERSE_PROFILE: {TYPE: UniformSphere, RMS_RADIUS: 0}\n");
  GaussianEmitter emitter(Flavour(kf_p_plus), 1.2, 1);
  emitter.SetERan(.5);
  CHECK(emitter.ImpactParameter() ==
        Catch::Approx(1.2 / rpa->hBar_c() * std::sqrt(std::log(2.))));
  CHECK(emitter.Weight() == 1.);
}

TEST_CASE("Actual EPA photon radii and weights ignore generic profile settings",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "EPA: {Form_Factor: Point-like, xBins: 4, bBins: 5, CacheTables: false}\n"
      "BEAM_TRANSVERSE_PROFILE:\n"
      "  TYPE: [Gaussian, UniformSphere]\n"
      "  RMS_RADIUS: [0, 99]\n");
  Settings::GetMainSettings()["GENERATE_RESULT_DIRECTORY"].SetDefault(false);
  Settings::GetMainSettings()["RESULT_DIRECTORY"].SetDefault("");
  EPA first(Flavour(kf_p_plus), 100., 0., 1);
  EPA second(Flavour(kf_p_plus), 100., 0., -1);
  for (double x : {1.e-4, .01, .2}) {
    for (double u : {.05, .4, .9}) {
      first.SetERan(u);
      second.SetERan(u);
      REQUIRE(first.CalculateWeight(x, 1.));
      REQUIRE(second.CalculateWeight(x, 1.));
      CHECK(first.ImpactParameter() == second.ImpactParameter());
      CHECK(first.Weight() == second.Weight());
      CHECK(std::isfinite(first.Weight()));
      CHECK(first.Weight() >= 0.);
    }
  }
}

TEST_CASE("Active overlap geometry with unit survival has unit weight",
          "[BEAM::Transverse_Profile]")
{
  for (const std::string& type : {"Gaussian", "UniformSphere", "WoodsSaxon"}) {
    BeamSettings settings(
        "BEAM_OVERLAP_REJECTION: 1\nBEAM_TRANSVERSE_PROFILE:\n"
        "  TYPE: " + type + "\n  RMS_RADIUS: 2.\n");
    // Electron beam radii are zero: the real radius rejector is identically
    // one, even though both transverse sampling distributions are extended.
    GaussianEmitter emitter(Flavour(kf_e), 1.1, 1);
    Monochromatic target(Flavour(kf_e), 100., 0., -1);
    TestKinematics kinematics(&emitter, &target);
    Collider_Weight weight(&kinematics);
    for (double u : {0., .2, .8, .999}) {
      emitter.SetERan(u);
      target.SetERan(1. - .5 * (1. + u));
      weight.SetERan(.37);
      CHECK(weight.OverlapWeight() == 1.);
      CHECK(weight() == 1.); // no second density, Jacobian, area or A factor
    }
  }
}

TEST_CASE("Gaussian sampling and actual overlap weight reproduce the analytic black disk",
          "[BEAM::Transverse_Profile]")
{
  const double r1_fm = 1.1, r2_fm = .9;
  BeamSettings settings(
      "BEAM_OVERLAP_REJECTION: 1\n"
      "BEAM_TRANSVERSE_PROFILE: {TYPE: Gaussian, RMS_RADIUS: 0.9}\n");
  GaussianEmitter emitter(Flavour(kf_p_plus), r1_fm, 1);
  Monochromatic target(Flavour(kf_p_plus), 100., 0., -1);
  TestKinematics kinematics(&emitter, &target);
  Collider_Weight weight(&kinematics);
  Monochromatic swapped_target(Flavour(kf_p_plus), 100., 0., 1);
  GaussianEmitter swapped_emitter(Flavour(kf_p_plus), r1_fm, -1);
  TestKinematics swapped_kinematics(&swapped_target, &swapped_emitter);
  Collider_Weight swapped_weight(&swapped_kinematics);

  const unsigned int n = 131072;
  double sum = 0., swap_difference = 0.;
  for (unsigned int i = 1; i <= n; ++i) {
    const double u = RadicalInverse(i, 2), v = RadicalInverse(i, 3);
    const double phi = RadicalInverse(i, 5);
    emitter.SetERan(u);
    target.SetERan(v);
    weight.SetERan(phi);
    const double w = weight();
    sum += w;
    if (i <= 256) {
      swapped_emitter.SetERan(u);
      swapped_target.SetERan(v);
      swapped_weight.SetERan(phi);
      swap_difference = std::max(swap_difference, std::abs(w - swapped_weight()));
    }
  }
  const double cut_fm = 2. * Flavour(kf_p_plus).Radius();
  const double expected = std::exp(-cut_fm * cut_fm /
                                   (r1_fm * r1_fm + r2_fm * r2_fm));
  CHECK(sum / n == Catch::Approx(expected).margin(8.e-4));
  CHECK(swap_difference == 0.);
}

TEST_CASE("Woods-Saxon projected CDF agrees with independent quadrature",
          "[BEAM::Transverse_Profile]")
{
  // SciPy quad/brentq references for rho(r)=1/(1+exp((r-R0)/a)), without
  // the production core/tail split, radius solve or interpolation table.
  struct Reference {
    double rms, skin;
    std::vector<std::pair<double, double>> quantiles;
  };
  const Reference references[] = {
      {4.49856217882625, .54,
       {{1.e-6, .00558303047239511}, {.01, .558813293415187},
        {.1, 1.78235435333479}, {.5, 4.19052915471790},
        {.9, 6.38138380776849}, {.99, 7.96185754963272},
        {.9999, 10.6849746183198}, {.999999, 13.3420095719641}}},
      {1., .3,
       {{1.e-6, .000929077310520787}, {.01, .0932176610673695},
        {.1, .302764152311552}, {.5, .787082418393334},
        {.9, 1.52823364991961}, {.99, 2.39102466382323},
        {.9999, 3.97626791173109}, {.999999, 5.49218044643381}}}};
  for (const auto& ref : references) {
    Transverse_Profile profile(Type::WoodsSaxon, ref.rms, ref.skin);
    CHECK(profile.Type() == Type::WoodsSaxon);
    CHECK(profile.SkinDepth() == ref.skin);
    CHECK(profile.SampleRadius(0.) == 0.);
    for (const auto& point : ref.quantiles) {
      const double u = point.first, radius = point.second;
      INFO("RMS=" << ref.rms << ", skin=" << ref.skin << ", u=" << u);
      // Bound the CDF error, not the radius error: the inverse is steep in
      // the dilute tail, so a uniform radius tolerance would be misleading.
      CHECK(profile.SampleRadius(u - 3.e-8) <= radius);
      CHECK(profile.SampleRadius(u + 3.e-8) >= radius);
    }
    const double edge = profile.SampleRadius(std::nextafter(1., 0.));
    CHECK(std::isfinite(edge));
    CHECK(edge > ref.quantiles.back().second);
  }
}

TEST_CASE("Woods-Saxon sampling reproduces the requested RMS and fourth moment",
          "[BEAM::Transverse_Profile]")
{
  const double rms = 4.49856217882625;
  Transverse_Profile profile(Type::WoodsSaxon, rms, .54);
  const unsigned int n = 200000;
  double second = 0., fourth = 0., previous = 0.;
  bool monotonic = true;
  for (unsigned int i = 0; i < n; ++i) {
    const double s = profile.SampleRadius((i + .5) / n);
    monotonic = monotonic && s >= previous;
    previous = s;
    second += s * s;
    fourth += s * s * s * s;
  }
  CHECK(monotonic);
  CHECK(second / n == Catch::Approx(rms * rms).epsilon(1.e-5));
  CHECK(fourth / n == Catch::Approx(632.795629221091).epsilon(1.e-4));

  // An ordinary value copy owns an equivalent table, with no shared state.
  const auto copy = profile;
  CHECK(copy.SampleRadius(.73) == profile.SampleRadius(.73));
}

TEST_CASE("Woods-Saxon has the pointlike and uniform-sphere limits",
          "[BEAM::Transverse_Profile]")
{
  Transverse_Profile point(Type::WoodsSaxon, 0., .54);
  Transverse_Profile sharp(Type::WoodsSaxon, 1., 0.);
  Transverse_Profile thin(Type::WoodsSaxon, 1., 1.e-4);
  Transverse_Profile sphere(Type::UniformSphere, 1.);
  for (double u : {0., .1, .5, .9, std::nextafter(1., 0.)}) {
    CHECK(point.SampleRadius(u) == 0.);
    CHECK(sharp.SampleRadius(u) == sphere.SampleRadius(u));
  }
  for (double u : {.1, .5, .9})
    CHECK(thin.SampleRadius(u) == Catch::Approx(sphere.SampleRadius(u)).epsilon(3.e-6));
}

TEST_CASE("Woods-Saxon rejects invalid or incompatible diffuseness",
          "[BEAM::Transverse_Profile]")
{
  Test_Setup::InitGlobals();
  for (double bad : {-1., std::numeric_limits<double>::infinity(),
                     std::numeric_limits<double>::quiet_NaN()})
    CHECK_THROWS(Transverse_Profile(Type::WoodsSaxon, 4.5, bad));
  CHECK_THROWS(Transverse_Profile(Type::WoodsSaxon, 1., .35));
  CHECK_THROWS(Transverse_Profile(Type::WoodsSaxon, 1., 1.));
  CHECK_THROWS(Transverse_Profile(Type::WoodsSaxon, 1.e-300, 1.e300));
  // R0=0 is a valid limiting profile: a/RMS=sqrt(zeta(3)/(10*zeta(5))).
  CHECK_NOTHROW(Transverse_Profile(Type::WoodsSaxon, 1., .34047737955445977478));
}

TEST_CASE("Woods-Saxon beam settings keep the common RMS convention and units",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_TRANSVERSE_PROFILE:\n"
      "  TYPE: [Gaussian, WoodsSaxon]\n"
      "  RMS_RADIUS: [0, 4.5]\n");
  Monochromatic target(Flavour(kf_p_plus), 100., 0., -1);
  CHECK(target.TransverseProfile().Type() == Type::WoodsSaxon);
  CHECK(target.TransverseProfile().RMSRadius() == Catch::Approx(4.5 / rpa->hBar_c()));
  CHECK(target.TransverseProfile().SkinDepth() == Catch::Approx(.54 / rpa->hBar_c()));
  Transverse_Profile reference(Type::WoodsSaxon, 4.5, .54);
  target.SetERan(.8);
  CHECK(target.ImpactParameter() * rpa->hBar_c() ==
        Catch::Approx(reference.SampleRadius(.8)).epsilon(1.e-7));
  std::unique_ptr<Beam_Base> copy(target.Copy());
  copy->SetERan(.8);
  CHECK(copy->ImpactParameter() == target.ImpactParameter());
}

TEST_CASE("Woods-Saxon skin depth broadcasts and accepts per-beam overrides",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_TRANSVERSE_PROFILE: {TYPE: WoodsSaxon, RMS_RADIUS: 4.5}\n");
  SECTION("per-beam overrides") {
    Settings::GetMainSettings()["BEAM_TRANSVERSE_PROFILE"]["SKIN_DEPTH"]
        .OverrideVector(std::vector<double>{.4, .6});
    Monochromatic first(Flavour(kf_p_plus), 100., 0., 1);
    Monochromatic second(Flavour(kf_p_plus), 100., 0., -1);
    CHECK(first.TransverseProfile().SkinDepth() == Catch::Approx(.4 / rpa->hBar_c()));
    CHECK(second.TransverseProfile().SkinDepth() == Catch::Approx(.6 / rpa->hBar_c()));
  }
  SECTION("scalar zero gives the analytic sphere limit on both beams") {
    Settings::GetMainSettings()["BEAM_TRANSVERSE_PROFILE"]["SKIN_DEPTH"]
        .OverrideVector(std::vector<double>{0.});
    const Transverse_Profile sphere(Type::UniformSphere, 4.5 / rpa->hBar_c());
    for (int dir : {1, -1}) {
      Monochromatic beam(Flavour(kf_p_plus), 100., 0., dir);
      beam.SetERan(.6);
      CHECK(beam.TransverseProfile().SkinDepth() == 0.);
      CHECK(beam.ImpactParameter() == sphere.SampleRadius(.6));
    }
  }
}

TEST_CASE("Woods-Saxon overlap agrees with independent spatial integration",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings(
      "BEAM_OVERLAP_REJECTION: 1\n"
      "BEAM_TRANSVERSE_PROFILE: {TYPE: WoodsSaxon, RMS_RADIUS: 4.491711}\n");
  if (s_kftable.find(kf_lead208) == s_kftable.end())
    AddParticle(kf_lead208, 193.75, 5.5012, 0., 246, 0, true, 1, "Pb208", "Pb208");
  GaussianEmitter emitter(Flavour(kf_lead208), 8., 1);
  Monochromatic target(Flavour(kf_lead208), 100., 0., -1);
  TestKinematics kinematics(&emitter, &target);
  Collider_Weight weight(&kinematics);
  // Hold b_gamma=8 fm fixed. The independent reference integrates the
  // longitudinal WS density and the exact angular survival fraction at
  // Bcut=2*5.5012 fm, for R0=6.61236651205933 fm, a=0.54 fm.
  emitter.SetERan(-std::expm1(-1.));
  const unsigned int n = 131072;
  double sum = 0.;
  for (unsigned int i = 1; i <= n; ++i) {
    target.SetERan(RadicalInverse(i, 2));
    weight.SetERan(RadicalInverse(i, 3));
    sum += weight();
  }
  CHECK(sum / n == Catch::Approx(.246528499636890).margin(8.e-4));
}

TEST_CASE("Ion radii are read as uniform spheres for the default width",
          "[BEAM::Transverse_Profile]")
{
  BeamSettings settings("");
  if (s_kftable.find(kf_lead208) == s_kftable.end())
    AddParticle(kf_lead208, 193.75, 5.5012, 0., 246, 0, true, 1, "Pb208", "Pb208");
  Monochromatic lead(Flavour(kf_lead208), 100., 0., 1);
  CHECK(lead.TransverseProfile().RMSRadius() ==
        Catch::Approx(std::sqrt(2. / 5.) * lead.Beam().Radius() / rpa->hBar_c()));
}

TEST_CASE("Optical Glauber overlap agrees with independent integration",
          "[BEAM::Beam_Rejection]")
{
  BeamSettings settings("");
  if (s_kftable.find(kf_lead208) == s_kftable.end())
    AddParticle(kf_lead208, 193.75, 5.5012, 0., 246, 0, true, 1, "Pb208", "Pb208");
  const Flavour pb(kf_lead208);
  const double hbarc = rpa->hBar_c(), sNN = sqr(5020.);
  Nucleus_Nucleus_Rejection rejector(pb, pb, 6.62 / hbarc, .546 / hbarc, sNN);

  // T_AA in 1/fm^2 from brute-force polar integration of the same density
  const std::vector<std::pair<double, double>> reference = {
      {0., 305.0015},    {5., 181.3027},      {10., 33.31511},
      {14., 1.016553},   {15., .2720558},     {15.7, .1006603},
      {17., .01431648},  {20., 1.195729e-4},  {23., 8.186443e-7}};
  const double sigmaNN = (25.0 + .146 * sqr(std::log(sNN))) * .1; // in fm^2
  for (const auto& [b, TAA] : reference) {
    CHECK(rejector.TAA(b / hbarc) / sqr(hbarc) ==
          Catch::Approx(TAA).epsilon(1.e-3));
    CHECK(rejector(b / hbarc) ==
          Catch::Approx(std::exp(-sigmaNN * TAA)).margin(5.e-4));
  }

  // the integral of T_AA over the transverse plane is A^2, and P(b) is monotonic
  const double db = .01 / hbarc;
  double integral = 0., last = 0.;
  for (double b = .5 * db; b < 60. / hbarc; b += db) {
    integral += 2. * M_PI * b * db * rejector.TAA(b);
    CHECK(rejector(b) >= last);
    last = rejector(b);
  }
  CHECK(integral == Catch::Approx(208. * 208.).epsilon(2.e-4));
  CHECK(last == 1.);
}

TEST_CASE("Identical nuclei use the optical Glauber overlap as event weight",
          "[BEAM::Beam_Rejection]")
{
  BeamSettings settings(
      "BEAM_OVERLAP_REJECTION: 2\n"
      "EPA: {Form_Factor: 16, xBins: 4, bBins: 5, CacheTables: false}\n"
      "BEAM_TRANSVERSE_PROFILE: {RMS_RADIUS: 4.5}\n");
  Settings::GetMainSettings()["GENERATE_RESULT_DIRECTORY"].SetDefault(false);
  Settings::GetMainSettings()["RESULT_DIRECTORY"].SetDefault("");
  if (s_kftable.find(kf_lead208) == s_kftable.end())
    AddParticle(kf_lead208, 193.75, 5.5012, 0., 246, 0, true, 1, "Pb208", "Pb208");
  const Flavour pb(kf_lead208);
  EPA emitter(pb, 100., 0., 1);
  Monochromatic target(pb, 100., 0., -1);
  TestKinematics kinematics(&emitter, &target);
  Collider_Weight weight(&kinematics);

  // the default density of the Woods-Saxon flux, whichever form factor is used
  const double hbarc = rpa->hBar_c();
  Nucleus_Nucleus_Rejection rejector(
      pb, pb, 1.118 * std::cbrt(208.) / hbarc, .54 / hbarc,
      (emitter.InMomentum() + target.InMomentum()).Abs2());
  bool transition = false;
  for (int i = 1; i < 100; ++i) {
    emitter.SetERan(.01 * i);
    REQUIRE(emitter.CalculateWeight(1.e-3, 1.));
    target.SetERan(.3);
    weight.SetERan(.2);
    const double b1 = emitter.ImpactParameter(), b2 = target.ImpactParameter();
    const double P = rejector(std::sqrt(
        b1 * b1 + b2 * b2 - 2. * b1 * b2 * std::cos(2. * M_PI * .2)));
    CHECK(weight.OverlapWeight() == Catch::Approx(P));
    CHECK(weight() == Catch::Approx(emitter.Weight() * target.Weight() * P));
    transition = transition || (P > .01 && P < .99);
  }
  CHECK(transition);

  // the other beam combinations fail at set-up
  if (s_kftable.find(kf_gold197) == s_kftable.end())
    AddParticle(kf_gold197, 183.5, 6.9823, 0., 237, 0, true, 1, "Au197", "Au197");
  GaussianEmitter gold(Flavour(kf_gold197), 15., 1);
  TestKinematics different(&gold, &target);
  CHECK_THROWS(Collider_Weight(&different));
  GaussianEmitter proton(Flavour(kf_p_plus), 15., 1);
  TestKinematics pA(&proton, &target);
  CHECK_THROWS(Collider_Weight(&pA));
}
