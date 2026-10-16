#include <catch2/catch_all.hpp>

#include "PHASIC++/Process/Massive_Kernels.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Phys/NLO_Types.H"

#include <cmath>
#include <vector>

using namespace ATOOLS;
using PHASIC::Massive_Kernels;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

  void Register(kf_code kf, double mass, int icharge, const std::string& name) {
    if (s_kftable.find(kf) == s_kftable.end())
      AddParticle(kf, mass, 0., 0., icharge, 1, true, 1, name, name);
  }

  void Boot() {
    static bool booted = false;
    if (booted) return;
    if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
    if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
    Settings::InitializeMainSettings("");
    // the same records as the other unit tests, so the shared table does not
    // depend on the test order; the q <- g kernels need neither colour nor masses
    Register(kf_d, 0.01, -1, "d");
    Register(kf_u, 0.005, 2, "u");
    Register(kf_s, 0.2, -1, "s");
    Register(kf_c, 1.42, 2, "c");
    Register(kf_b, 4.92, -1, "b");
    booted = true;
  }

  // SetNC sums Q^2 over the model's massless charged flavours; only the
  // quarks matter, and only for the I-operator constants, not the KP kernels.
  class Quark_Model : public MODEL::Model_Base {
  public:
    Quark_Model() : Model_Base(true) {
      for (long kf = 1; kf <= 5; ++kf) {
        m_fls.push_back(Flavour(kf));
        m_fls.push_back(Flavour(kf).Bar());
      }
    }
    void InitVertices() override {}
    void ParticleInit() override {}
    bool ModelInit() override { return true; }
  };

  // Massive_Kernels with N_C = 3 set the way KP_Terms does it
  struct Kernels {
    Quark_Model model;
    MODEL::Model_Base *saved;
    Massive_Kernels kernels;
    explicit Kernels(sbt::subtype st) : saved(MODEL::s_model), kernels(st, 5, 0) {
      MODEL::s_model = &model;
      kernels.SetNC(3.);
    }
    ~Kernels() { MODEL::s_model = saved; }
  };

  // regular part of the q <- g (QCD) or q <- photon (QED) splitting function
  // per T^2 of the quark entering the Born, without its colour/charge factor
  double Pqg(double x) { return x * x + (1. - x) * (1. - x); }

  const std::vector<double> xs = {0.05, 0.3, 0.7, 0.95};
  constexpr int q_from_g = 3;  // Massive_Kernels type: quark entering, gluon/photon split
  constexpr int quark = 1;     // spin of a final-state quark spectator

}

TEST_CASE("Quark-from-gluon and quark-from-photon KP kernels follow Catani-Seymour",
          "[PHASIC::Massive_Kernels]") {
  Boot();
  // T_R/C_F for g -> q qbar; for gamma -> q qbar T_R -> N_C Q^2, C_F -> Q^2
  const auto [st, c] = GENERATE(table<sbt::subtype, double>(
      {{sbt::qcd, 0.5 / (4. / 3.)}, {sbt::qed, 3.}}));
  {
    Kernels k(st);
    auto &K = k.kernels;
    CAPTURE(st, c);

    SECTION("an off-diagonal splitting has no endpoint or plus-distribution pieces") {
      K.SetAlpha(0.1, 0.1, 0.1, 0.1);
      CHECK(K.Kb2(q_from_g) == 0.);
      CHECK(K.Kt2(q_from_g) == 0.);
      CHECK(K.P2(q_from_g) == 0.);
      CHECK(K.t2(q_from_g, quark, 0.) == 0.);
      for (double x : xs) {
        CAPTURE(x);
        CHECK(K.Kb1(q_from_g, x) == 0.);
        CHECK(K.Kb4(q_from_g, x) == 0.);
        CHECK(K.Kt1(q_from_g, x) == 0.);
        CHECK(K.Kt4(q_from_g, x) == 0.);
        CHECK(K.P1(q_from_g, x) == 0.);
        CHECK(K.P4(q_from_g, x) == 0.);
        CHECK(K.t1(q_from_g, quark, 0., x) == 0.);
        CHECK(K.t4(q_from_g, quark, 0., x) == 0.);
      }
    }

    SECTION("Kbar is c [P(x) ln((1-x)/x) + 2x(1-x)] and P is c P(x)") {
      for (double x : xs) {
        CAPTURE(x);
        CHECK_THAT(K.Kb3(q_from_g, x),
                   WithinRel(c * (Pqg(x) * std::log((1. - x) / x) + 2. * x * (1. - x)), 1e-12));
        CHECK_THAT(K.P3(q_from_g, x), WithinRel(c * Pqg(x), 1e-12));
      }
    }

    SECTION("Ktilde is c P(x) ln(1-x) at alpha_II = 1") {
      for (double x : xs) {
        CAPTURE(x);
        CHECK_THAT(K.Kt3(q_from_g, x), WithinRel(c * Pqg(x) * std::log(1. - x), 1e-12));
      }
    }

    SECTION("alpha_II < 1 cuts Ktilde to c P(x) ln(alpha_II) where 1-x > alpha_II") {
      const double a = 0.1;
      K.SetAlpha(1., 1., 1., a);
      CHECK_THAT(K.Kt3(q_from_g, 0.3), WithinRel(c * Pqg(0.3) * std::log(a), 1e-12));
      CHECK_THAT(K.Kt3(q_from_g, 0.95),
                 WithinRel(c * Pqg(0.95) * std::log(1. - 0.95), 1e-12));
    }

    SECTION("in the CSS scheme Ktilde gains the integral of the II kernel shift, -c (1-x)(1-3x) at alpha_II = 1") {
      // the CSS II dipole evaluates 1-2z(1-z) at z = x+v instead of x (FFV_DLC);
      // int_0^{1-x} dv/v [-2v(1-2x-v)] = -(1-x)(1-3x)
      K.SetSubType(subscheme::CSS);
      for (double x : xs) {
        CAPTURE(x);
        CHECK_THAT(K.Kt3(q_from_g, x),
                   WithinRel(c * (Pqg(x) * std::log(1. - x) - (1. - x) * (1. - 3. * x)), 1e-12));
      }
    }

    SECTION("a massless final-state spectator adds nothing at alpha_IF = 1") {
      K.SetAlpha(1., 0.1, 1., 1.);  // alpha_FI restricts final-state emitters only
      for (double x : xs) {
        CAPTURE(x);
        CHECK_THAT(K.t3(q_from_g, quark, 0., x), WithinAbs(0., 1e-15));
      }
    }

    SECTION("alpha_IF < 1 adds the spectator term P(x) ln(1/alpha_IF) in the same normalisation as P") {
      const double a = 0.1;
      K.SetAlpha(1., 1., a, 1.);
      for (double x : xs) {
        CAPTURE(x);
        CHECK_THAT(K.t3(q_from_g, quark, 0., x),
                   WithinRel(K.P3(q_from_g, x) * std::log(1. / a), 1e-12));
      }
    }
  }
}
