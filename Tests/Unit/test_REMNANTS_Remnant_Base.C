#include <catch2/catch_all.hpp>

#include "REMNANTS/Main/Remnant_Base.H"
#include "REMNANTS/Tools/Colour_Generator.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace ATOOLS;
using namespace REMNANTS;

namespace {

  void Boot() {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Settings::InitializeMainSettings("");
    if (s_kftable.find(kf_photon) == s_kftable.end()) {
      AddParticle(kf_photon, 0., 0., 0., 0, 2, true, 1, "P", "P");
      AddParticle(kf_d, 0.01, 0., 0., -1, 1, true, 1, "d", "d");
    }
    booted = true;
  }

  class Test_Beam : public BEAM::Beam_Base {
  public:
    Test_Beam(const Flavour& flav, const double energy)
        : Beam_Base(BEAM::beamspectrum::monochromatic, flav, energy, 0., 1) {}
    Beam_Base* Copy() override { return nullptr; }
    bool CalculateWeight(double, double) override { return true; }
  };

  class Test_Remnant : public Remnant_Base {
  public:
    using Remnant_Base::Remnant_Base;
    bool TestExtract(const Flavour& flav, const Vec4D& mom) override {
      return true;
    }
    bool FillBlob(Colour_Generator* colours, ParticleMomMap* ktmap = nullptr,
                  const bool& copy = true) override {
      return true;
    }
  };

}

TEST_CASE("Residual equals the beam-out momentum before extraction",
          "[REMNANTS::Remnant_Base]") {
  Boot();
  Test_Beam beam(Flavour(kf_photon), 10.);
  Test_Remnant remnant(Flavour(kf_photon), 0);
  remnant.SetBeam(&beam);
  remnant.Reset();

  beam.SetOutMomentum(Vec4D(4., 0., 0., 4.));
  CHECK(remnant.Residual() == Vec4D(4., 0., 0., 4.));
  CHECK(remnant.Residual() == remnant.IncomingMomentum());
}

TEST_CASE("Residual tracks extracted partons", "[REMNANTS::Remnant_Base]") {
  Boot();
  Test_Beam beam(Flavour(kf_photon), 10.);
  Test_Remnant remnant(Flavour(kf_photon), 0);
  remnant.SetBeam(&beam);
  remnant.Reset();
  beam.SetOutMomentum(Vec4D(4., 0., 0., 4.));

  Colour_Generator colours;
  Particle part1(-1, Flavour(kf_d), Vec4D(1., 0., 0., 1.));
  Particle part2(-1, Flavour(kf_d).Bar(), Vec4D(0.5, 0., 0., 0.5));

  CHECK(remnant.Extract(&part1, &colours));
  CHECK(remnant.Residual() == Vec4D(3., 0., 0., 3.));

  CHECK(remnant.Extract(&part2, &colours));
  CHECK(remnant.Residual() == Vec4D(2.5, 0., 0., 2.5));

  // Re-extracting the same particle must not be double-counted.
  CHECK(remnant.Extract(&part1, &colours));
  CHECK(remnant.Residual() == Vec4D(2.5, 0., 0., 2.5));
}

TEST_CASE("Residual picks up a new beam-out momentum without Reset",
          "[REMNANTS::Remnant_Base]") {
  // Regression: the residual must follow per-event rewrites of the beam-out
  // momentum (e.g. EPA photons) on demand, not via any cached scalar.
  Boot();
  Test_Beam beam(Flavour(kf_photon), 10.);
  Test_Remnant remnant(Flavour(kf_photon), 0);
  remnant.SetBeam(&beam);
  remnant.Reset();
  beam.SetOutMomentum(Vec4D(4., 0., 0., 4.));

  Colour_Generator colours;
  Particle part(-1, Flavour(kf_d), Vec4D(1., 0., 0., 1.));
  CHECK(remnant.Extract(&part, &colours));
  CHECK(remnant.Residual() == Vec4D(3., 0., 0., 3.));

  beam.SetOutMomentum(Vec4D(6., 0., 0., 6.));
  CHECK(remnant.Residual() == Vec4D(5., 0., 0., 5.));

  remnant.Reset();
  CHECK(remnant.Residual() == Vec4D(6., 0., 0., 6.));
}
