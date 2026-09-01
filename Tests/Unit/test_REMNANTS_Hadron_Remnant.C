#include <catch2/catch_all.hpp>

#include "REMNANTS/Main/Hadron_Remnant.H"
#include "REMNANTS/Tools/Remnants_Parameters.H"
#include "PDF/Main/PDF_Base.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace ATOOLS;
using namespace REMNANTS;

namespace {

  void Register(kf_code kf, double mass, int icharge, int spin,
                const std::string& name) {
    if (s_kftable.find(kf) == s_kftable.end())
      AddParticle(kf, mass, 0., 0., icharge, spin, true, 1, name, name);
  }

  void Boot() {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Settings::InitializeMainSettings("");
    Register(kf_photon, 0., 0, 2, "P");
    Register(kf_gluon, 0., 0, 2, "G");
    Register(kf_d, 0.01, -1, 1, "d");
    Register(kf_u, 0.005, 2, 1, "u");
    Register(kf_e, 0.000511, -3, 1, "e-");
    Register(kf_mu, 0.1057, -3, 1, "mu-");
    Register(kf_pi, 0.135, 0, 0, "pi");
    Register(kf_pi_plus, 0.1396, 3, 0, "pi+");
    Register(kf_n, 0.9396, 0, 1, "n");
    Register(kf_p_plus, 0.938, 3, 1, "P+");
    if (!rempars) {
      rempars = new Remnants_Parameters();
      rempars->Init();
    }
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

TEST_CASE("Extraction checks are invariant under longitudinal boosts",
          "[REMNANTS::Hadron_Remnant]") {
  // HERA-like setup: 920 GeV proton along -z, low-energy photon partner
  // along +z.  The accept/reject decision must be the same whether the
  // momenta are given in the lab or in the photon-proton c.m. frame.
  Boot();
  Test_Beam beam(Flavour(kf_p_plus), 920., -1);
  Test_PDF pdf(Flavour(kf_p_plus));
  Hadron_Remnant remnant(&pdf, 1);
  remnant.SetBeam(&beam);
  remnant.Reset();

  const Vec4D pphoton(4., 0., 0., 4.);
  const Vec4D pproton = beam.OutMomentum();
  const double spair = (pphoton + pproton).Abs2();
  Poincare cms(pphoton + pproton);

  const std::vector<Vec4D> candidates = {
    Vec4D(300., 1., 0., -299.9),   // moderate x, plenty of energy left
    Vec4D(919., 0.5, 0., -918.9),  // takes almost the full proton
    Vec4D(0.001, 0., 0., -0.001),  // below the PDF x range
  };
  for (const Vec4D& mom : candidates) {
    const bool lab = remnant.TestExtract(Flavour(kf_u), mom, spair);
    beam.SetOutMomentum(cms * pproton);
    const bool boosted = remnant.TestExtract(Flavour(kf_u), cms * mom, spair);
    beam.SetOutMomentum(pproton);
    CHECK(lab == boosted);
  }
}

TEST_CASE("Energy reserve is evaluated in the pair c.m. frame",
          "[REMNANTS::Hadron_Remnant]") {
  Boot();
  Test_Beam beam(Flavour(kf_p_plus), 920., -1);
  Test_PDF pdf(Flavour(kf_p_plus));
  Hadron_Remnant remnant(&pdf, 1);
  remnant.SetBeam(&beam);
  remnant.Reset();

  // sqrt(spair) ~ 121 GeV for a 4 GeV photon on a 920 GeV proton.
  const Vec4D pphoton(4., 0., 0., 4.);
  const double spair = (pphoton + beam.OutMomentum()).Abs2();

  // x ~ 0.5: the leftover carries ~ sqrt(spair)/4 ~ 30 GeV in the pair
  // c.m. frame, well above the reserve of m_minE = 2 GeV for baryons.
  CHECK(remnant.TestExtract(Flavour(kf_u), Vec4D(460., 0., 0., -460.), spair));

  // x ~ 0.99: the leftover carries ~ 9 GeV in the lab, but only
  // ~ 0.01*sqrt(spair)/2 ~ 0.6 GeV in the pair c.m. frame - rejected.
  CHECK(!remnant.TestExtract(Flavour(kf_u), Vec4D(910., 0., 0., -910.), spair));
}
