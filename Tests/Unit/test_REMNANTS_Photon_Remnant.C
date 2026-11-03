#include <catch2/catch_all.hpp>

#include "Test_Setup.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "PDF/Main/PDF_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Org/Settings.H"

using namespace ATOOLS;
using namespace REMNANTS;

namespace {

  void Boot() {
    static bool booted = false;
    if (booted) return;
    Test_Setup::InitGlobals();
    Settings::InitializeMainSettings("");
    Test_Setup::RegisterParticles();
    if (!rempars) rempars = new Remnants_Parameters();
    booted = true;
  }

  class Test_Beam : public BEAM::Beam_Base {
  public:
    Test_Beam(const Flavour& flav, const double energy, const int dir)
        : Beam_Base(BEAM::beamspectrum::monochromatic, flav, energy, 0., dir) {}
    Beam_Base* Copy() override { return nullptr; }
    bool CalculateWeight(double, double) override { return true; }
  };

  class Test_PDF : public PDF::PDF_Base {
  public:
    explicit Test_PDF(const Flavour& bunch) {
      m_bunch = bunch;
      m_partons.insert(Flavour(kf_d));
      m_partons.insert(Flavour(kf_d).Bar());
      m_partons.insert(Flavour(kf_u));
      m_partons.insert(Flavour(kf_u).Bar());
      m_partons.insert(Flavour(kf_gluon));
      m_xmin  = 1.e-5;
      m_xmax  = 1.;
      m_q2min = 1.;
      m_q2max = 1.e6;
      SetBounds();
    }
    void CalculateSpec(const double&, const double&) override {}
    double GetXPDF(const ATOOLS::Flavour&) override { return 1.; }
    double GetXPDF(const kf_code&, bool) override { return 1.; }
    PDF_Base* GetCopy() override { return nullptr; }
  };

}

TEST_CASE("Momentum fraction uses light-cone components",
          "[REMNANTS::Photon_Remnant]") {
  // 4 GeV EPA-like photon along +z.  The extraction fraction is defined on
  // light-cone components, so a candidate with large transverse momentum can
  // have an energy above the photon energy and must still be accepted, while
  // the old energy fraction would have put it beyond the PDF x range.
  Boot();
  Test_Beam beam(Flavour(kf_photon), 4., 1);
  Test_PDF pdf(Flavour(kf_photon));
  Photon_Remnant remnant(&pdf, 0);
  remnant.SetBeam(&beam);
  remnant.Reset();

  const double spair = (beam.OutMomentum() + Vec4D(920., 0., 0., -920.)).Abs2();

  // E/E_photon = 1.125, but x_lc = PPlus/PPlus_photon = 4.9/8 ~ 0.61.
  CHECK(remnant.TestExtract(Flavour(kf_u), Vec4D(4.5, 4.4, 0., 0.4), spair));

  // A backward-moving candidate has no forward light-cone momentum.
  CHECK(!remnant.TestExtract(Flavour(kf_u), Vec4D(1., 0., 0., -1.), spair));

  // Below the PDF x range.
  CHECK(!remnant.TestExtract(Flavour(kf_u), Vec4D(1.e-6, 0., 0., 1.e-6), spair));
}
