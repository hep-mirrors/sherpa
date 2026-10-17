#include <catch2/catch_all.hpp>

#include "PHASIC++/Process/KP_Terms.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PDF/Main/PDF_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "PDF/Main/Structure_Function.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Coupling_Data.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Phys/Weights.H"

#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

using namespace ATOOLS;
using PHASIC::KP_Terms;
using Catch::Matchers::WithinRel;

CATCH_TRANSLATE_EXCEPTION(const ATOOLS::Exception& e) {
  std::ostringstream str;
  str << e;
  return str.str();
}

namespace {

  void Register(kf_code kf, double mass, int icharge, int spin,
                const std::string& name) {
    if (s_kftable.find(kf) == s_kftable.end())
      AddParticle(kf, mass, 0., 0., icharge, spin, true, 1, name, name);
  }

  // the settings KP_Terms reads, registered as Sherpa's initialisation and
  // the ME generators do, except for the subtraction scheme: plain
  // Catani-Seymour (Sherpa's default is CSS)
  void Seed(const std::string& yaml = "") {
    Settings::InitializeMainSettings(yaml);
    Settings& s = Settings::GetMainSettings();
    s["N_COLOR"].SetDefault(3.0);
    s["DIPOLES"]["SCHEME"].SetDefault(subscheme::CS);
    s["DIPOLES"]["KAPPA"].SetDefault(2.0 / 3.0);
    for (const char* a : {"ALPHA", "ALPHA_FF", "ALPHA_FI", "ALPHA_IF", "ALPHA_II"})
      s["DIPOLES"][a].SetDefault(1.0);
    s["DIPOLES"]["COLLINEAR_VFF_SPLITTINGS"].SetDefault(1);
    s["DIPOLES"]["V_SUBTRACTION_MODE"].SetDefault(1);
  }

  void Boot() {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Seed();
    // the same records as the REMNANTS tests, so the shared table does not
    // depend on the order the test files run in
    Register(kf_photon, 0., 0, 2, "P");
    Register(kf_gluon, 0., 0, 2, "G");
    Register(kf_d, 0.01, -1, 1, "d");
    Register(kf_u, 0.005, 2, 1, "u");
    Register(kf_p_plus, 0.938, 3, 1, "P+");
    // KP_Terms only counts the light flavours in this container; it cannot
    // hold both u and d here, as the shared records give them different masses
    if (s_kftable.find(kf_quark) == s_kftable.end()) {
      AddParticle(kf_quark, 0., 0., 0., 0, 1, 1, 0, 1, 1, 0, "Q", "Q", "Q", "Q", 1, 1);
      s_kftable[kf_quark]->Clear();
      s_kftable[kf_quark]->Add(Flavour(kf_u));
      s_kftable[kf_quark]->Add(Flavour(kf_u).Bar());
    }
    booted = true;
  }

  // Massive_Kernels::SetNC reads the model's flavours
  class Test_Model : public MODEL::Model_Base {
  public:
    Test_Model() : Model_Base(true) {
      for (kf_code kf : {kf_d, kf_u}) {
        m_fls.push_back(Flavour(kf));
        m_fls.push_back(Flavour(kf).Bar());
      }
    }
    void InitVertices() override {}
    void ParticleInit() override {}
    bool ModelInit() override { return true; }
  };

  class Fixed_Scale : public PHASIC::Scale_Setter_Base {
  public:
    explicit Fixed_Scale(double muf2) : Scale_Setter_Base(Args()) {
      m_scale[stp::fac] = muf2;
      m_fac.assign(stp::size, 1.);
    }
  protected:
    double Calculate(const Vec4D_Vector&, const size_t&) override
    { return m_scale[stp::fac]; }
  private:
    static PHASIC::Scale_Setter_Arguments Args() {
      PHASIC::Scale_Setter_Arguments args(nullptr);
      args.m_nin = 2;
      args.m_nout = 2;
      return args;
    }
  };

  class Test_Process : public PHASIC::Process_Base {
  public:
    Test_Process(const Flavour_Vector& flavs, double muf2) {
      m_flavs = flavs;
      m_nin = 2;
      m_nout = flavs.size() - 2;
      p_scale = new Fixed_Scale(muf2);
    }
    using Process_Base::Differential;
    void SetScale(const PHASIC::Scale_Setter_Arguments&) override {}
    void SetKFactor(const PHASIC::KFactor_Setter_Arguments&) override {}
    void InitializeTheReweighting(Variations_Mode) override {}
    size_t Size() const override { return 1; }
    Process_Base* operator[](const size_t&) override { return this; }
    Weight_Info* OneEvent(const int, Variations_Mode, const int) override
    { return nullptr; }
    Weights_Map Differential(const Vec4D_Vector&, Variations_Mode) override
    { return Weights_Map{}; }
    bool CalculateTotalXSec(const std::string&, const bool) override
    { return false; }
    void SetLookUp(const bool) override {}
  };

  // positive shapes without a photon; the pointlike term must not depend on them
  class Test_PDF : public PDF::PDF_Base {
    double m_x{0.5};
  public:
    explicit Test_PDF(const Flavour& bunch) {
      m_bunch = bunch;
      for (kf_code kf : {kf_d, kf_u}) {
        m_partons.insert(Flavour(kf));
        m_partons.insert(Flavour(kf).Bar());
      }
      m_partons.insert(Flavour(kf_gluon));
      m_xmin = 1.e-5;
      m_xmax = 1.;
      m_q2min = 1.;
      m_q2max = 1.e8;
      SetBounds();
    }
    void CalculateSpec(const double& x, const double&) override { m_x = x; }
    double GetXPDF(const Flavour& fl) override
    { return XPDF(fl.Kfcode(), fl.IsAnti()); }
    double GetXPDF(const kf_code& kf, bool anti) override { return XPDF(kf, anti); }
    PDF_Base* GetCopy() override { return nullptr; }
    double XPDF(kf_code kf, bool anti) const {
      const double n = kf == kf_gluon ? 1. : (anti ? 0.1 : 0.3) * double(kf);
      return n * std::pow(m_x, 0.3) * std::pow(1. - m_x, 2.);
    }
  };

  // massless 2->2 Born: parton 0 with energy e0 along dir*z, parton 1 with
  // energy e1 against it, partonic scattering angle theta
  Vec4D_Vector Born(double e0, double e1, double theta, int dir) {
    const double h = std::sqrt(e0 * e1), y = 0.5 * std::log(e0 / e1) * dir;
    auto lab = [&](double px, double pz) {
      return Vec4D(h * std::cosh(y) + pz * std::sinh(y), px, 0.,
                   pz * std::cosh(y) + h * std::sinh(y));
    };
    return {Vec4D(e0, 0., 0., dir * e0), Vec4D(e1, 0., 0., -dir * e1),
            lab(h * std::sin(theta), dir * h * std::cos(theta)),
            lab(-h * std::sin(theta), -dir * h * std::cos(theta))};
  }

  double Pqg(double x) { return x * x + (1. - x) * (1. - x); }

  // Catani-Seymour gamma -> q qbar pointlike term per unit Born, QED
  // replacement T_R -> N_C Q^2, T^2 -> Q^2, T.T -> charge correlators:
  //   alpha/2pi N_C Q_a^2 [ Kbar(eta) + w_b P ln(1-eta)
  //                         - P sum_k w_k ln(muF^2/s_0k)
  //                         - P ln(1/alpha_IF) sum_{k final} w_k ]
  // with spectator weights w_k = -eta_0 eta_k Q_0 Q_k / Q_0^2, sum_k w_k = 1
  struct Spectator { size_t index; double w; };
  double Pointlike(const Vec4D_Vector& p, double qa, double eta, double muf2,
                   const std::vector<Spectator>& spectators, double alpha_if = 1.) {
    double res = Pqg(eta) * std::log((1. - eta) / eta) + 2. * eta * (1. - eta);
    for (const auto& k : spectators) {
      res -= k.w * Pqg(eta) * std::log(muf2 / (2. * (p[0] * p[k.index])));
      if (k.index == 1) res += k.w * Pqg(eta) * std::log(1. - eta);
      else              res -= k.w * Pqg(eta) * std::log(1. / alpha_if);
    }
    const double alpha = 1. / 137.036;
    return alpha / (2. * M_PI) * 3. * qa * qa * res;
  }

  struct Setup {
    Test_Model model;
    MODEL::Model_Base* saved_model;
    Flavour saved_bunch[2];
    Test_Process proc;
    MODEL::Coupling_Data aqed;
    MODEL::Coupling_Map cpls;
    // beam-ordered; the ISR handler owns the PDFs
    Test_PDF *photon, *proton;
    std::unique_ptr<PDF::ISR_Handler> isr;
    Setup(const Flavour_Vector& flavs, double muf2, int dir)
        : saved_model(MODEL::s_model),
          saved_bunch{rpa->gen.Bunch(0), rpa->gen.Bunch(1)},
          proc(flavs, muf2), aqed(nullptr, "Alpha_QED", 1. / 137.036),
          photon(new Test_PDF(Flavour(kf_photon))),
          proton(new Test_PDF(Flavour(kf_p_plus))) {
      MODEL::s_model = &model;
      // the photon bunch sits on the +z side for dir=+1 and on -z for dir=-1
      Test_PDF* beam[2] = {dir > 0 ? photon : proton, dir > 0 ? proton : photon};
      isr.reset(new PDF::ISR_Handler(
          {new PDF::Structure_Function(beam[0], beam[0]->Bunch()),
           new PDF::Structure_Function(beam[1], beam[1]->Bunch())}));
      proc.Integrator()->SetISR(isr.get());
      for (size_t i = 0; i < 2; ++i) rpa->gen.SetBunch(beam[i]->Bunch(), i);
      cpls.insert(std::make_pair(std::string("Alpha_QED"), &aqed));
    }
    ~Setup() {
      MODEL::s_model = saved_model;
      rpa->gen.SetBunch(saved_bunch[0], 0);
      rpa->gen.SetBunch(saved_bunch[1], 1);
    }
    // Get() x f_a(eta0): the pointlike coefficient per unit Born for the
    // flavour-ordered momenta p, as the ME generators combine it with the
    // Born PDFs
    double Coefficient(const Vec4D_Vector& p, double eta0, double eta1,
                       double muf2, double alpha_if = 1.) {
      std::unique_ptr<KP_Terms> kp(
          KP_Terms::PhotonSplitting(&proc, &cpls, cs_itype::K | cs_itype::P));
      REQUIRE(kp);
      kp->SetAlpha(1., 1., alpha_if, 1.);
      const double x0 = 0.5 * (1. + eta0), x1 = 0.5 * (1. + eta1);
      const std::vector<std::vector<double>> dsij(1, std::vector<double>(1, 1.));
      kp->Calculate(p, dsij, x0, x1, eta0, eta1);
      // Get() is flavour-ordered: parton 0 comes from the photon
      const Flavour fl0 = proc.Flavours()[0], fl1 = proc.Flavours()[1];
      const double res = kp->Get(photon, proton, x0, x1, eta0, eta1,
                                 muf2, muf2, 1., 1., fl0, fl1);
      photon->Calculate(eta0, muf2);
      return res * photon->GetXPDF(fl0) / eta0;
    }
  };

  const double muf2 = 91.1876 * 91.1876;
  const double eta0 = 0.4, eta1 = 0.05;

  // seeds the settings with extra YAML and restores the defaults on exit
  struct Seeded {
    explicit Seeded(const std::string& yaml) { Seed(yaml); }
    ~Seeded() { Seed(); }
  };

}

TEST_CASE("The resolved-photon pointlike KP term is the Catani-Seymour gamma -> q qbar term",
          "[PHASIC::KP_Terms]") {
  Boot();
  const Flavour u(kf_u), d(kf_d), g(kf_gluon);
  const int dir = GENERATE(1, -1);
  CAPTURE(dir);

  SECTION("u g -> u g: the final-state u is the only spectator, with weight 1") {
    Setup s({u, g, u, g}, muf2, dir);
    const Vec4D_Vector p = Born(10., 41., 0.7, dir);
    CHECK_THAT(s.Coefficient(p, eta0, eta1, muf2),
               WithinRel(Pointlike(p, 2. / 3., eta0, muf2, {{2, 1.}}), 1e-6));
  }

  SECTION("u dbar -> u dbar: spectator weights -1/2 (initial dbar), 1 (u), 1/2 (dbar); the proton side adds nothing") {
    Setup s({u, d.Bar(), u, d.Bar()}, muf2, dir);
    const Vec4D_Vector p = Born(10., 41., 1.1, dir);
    CHECK_THAT(s.Coefficient(p, eta0, eta1, muf2),
               WithinRel(Pointlike(p, 2. / 3., eta0, muf2,
                                   {{1, -0.5}, {2, 1.}, {3, 0.5}}), 1e-6));
  }

  SECTION("alpha_IF < 1 adds -P ln(1/alpha_IF) times the final-state spectator weights") {
    Setup s({u, d.Bar(), u, d.Bar()}, muf2, dir);
    const Vec4D_Vector p = Born(10., 41., 1.1, dir);
    const double a = 0.1;
    CHECK_THAT(s.Coefficient(p, eta0, eta1, muf2, a),
               WithinRel(Pointlike(p, 2. / 3., eta0, muf2,
                                   {{1, -0.5}, {2, 1.}, {3, 0.5}}, a), 1e-6));
  }
}

TEST_CASE("In the DISgamma scheme no pointlike KP term is built", "[PHASIC::KP_Terms]") {
  Boot();
  const Seeded disgamma("KP: {FACTORISATION_SCHEME: DISgamma}\n");
  const Flavour u(kf_u), g(kf_gluon);
  Setup s({u, g, u, g}, muf2, 1);
  std::unique_ptr<KP_Terms> kp(
      KP_Terms::PhotonSplitting(&s.proc, &s.cpls, cs_itype::K | cs_itype::P));
  CHECK_FALSE(kp);
}
