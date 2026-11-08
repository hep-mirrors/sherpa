#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "YFS/Main/YFS_Form_Factor.H"

#include <cmath>
#include <memory>
#include <vector>

using ATOOLS::Vec4D;
using Catch::Approx;

namespace {

  const double alpha0 = 1./137.03599976;

  // YFS_Base reads alpha_QED from the model; nothing else of it is needed
  class Alpha_Only_Model : public MODEL::Model_Base {
  public:
    Alpha_Only_Model() : Model_Base(true) { (*p_constants)["alpha_QED"] = alpha0; }
    void InitVertices() override {}
    void ParticleInit() override {}
    bool ModelInit() override { return true; }
  };

  YFS::YFS_Form_Factor &FormFactor()
  {
    static std::unique_ptr<YFS::YFS_Form_Factor> ff;
    if (!ff) {
      if (!ATOOLS::msg) ATOOLS::msg = new ATOOLS::Message();
      if (!ATOOLS::rpa) ATOOLS::rpa = new ATOOLS::Run_Parameter();
      ATOOLS::Settings::InitializeMainSettings("");
      if (!MODEL::s_model) MODEL::s_model = new Alpha_Only_Model();
      ff.reset(new YFS::YFS_Form_Factor());
    }
    return *ff;
  }

  Vec4D OnShell(double m, double px, double py, double pz)
  {
    return Vec4D(std::sqrt(m*m + px*px + py*py + pz*pz), px, py, pz);
  }

  struct Pair {
    const char *name;
    Vec4D p1, p2;
    double m1, m2;
  };

  // back to back in the pair's rest frame, then boosted along a common vector
  Pair MakePair(const char *name, double sqrts, double m1, double m2,
                const Vec4D &boost = Vec4D(1., 0., 0., 0.))
  {
    const double lam(std::sqrt((sqrts*sqrts - (m1 + m2)*(m1 + m2))
                               *(sqrts*sqrts - (m1 - m2)*(m1 - m2))));
    const double p(lam/(2.*sqrts));
    Vec4D p1(OnShell(m1, 0.3*p, 0., std::sqrt(0.91)*p));
    Vec4D p2(OnShell(m2, -0.3*p, 0., -std::sqrt(0.91)*p));
    ATOOLS::Poincare L(boost);
    L.BoostBack(p1);
    L.BoostBack(p2);
    return Pair{name, p1, p2, m1, m2};
  }

  std::vector<Pair> Pairs()
  {
    return {
      MakePair("e+e- at the Z pole", 91.1876, 0.000510999, 0.000510999),
      MakePair("mu+mu- at 10.58 GeV, boosted", 10.58, 0.105658, 0.105658,
               Vec4D(12., 1., 2., 5.)),
      MakePair("pi+pi- at 1 GeV", 1., 0.13957, 0.13957),
      MakePair("e mu, unequal masses, boosted", 30., 0.000510999, 0.105658,
               Vec4D(40., -3., 0., 20.)),
    };
  }

  double Real(const Pair &x, double kmax, double lambda)
  {
    return FormFactor().BVR_full(x.p1*x.p2, x.p1[0], x.p2[0], x.m1, x.m2,
                                 kmax, lambda, 0);
  }

  // explicit masses: p.Mass() of a 45 GeV electron has lost six digits
  double Virtual(const Pair &x, double lambda)
  {
    return FormFactor().BVV_full_impl(x.p1, x.p2, x.m1, x.m2, lambda, 0);
  }

}

TEST_CASE("YFS form factor: real + virtual does not depend on the photon mass",
          "[YFS][FormFactor]") {
  for (const Pair &x : Pairs()) {
    INFO(x.name);
    const double kmax(0.01*(x.p1[0] + x.p2[0]));
    const double ref(Real(x, kmax, 1e-10) + Virtual(x, 1e-10));
    for (double lambda : {1e-20, 1e-15, 1e-6, 1e-4}) {
      INFO("lambda = " << lambda);
      const double r(Real(x, kmax, lambda)), v(Virtual(x, lambda));
      // each piece does depend on it
      CHECK(std::abs(r - Real(x, kmax, 1e-10)) > 1e-3);
      CHECK(r + v == Approx(ref).margin(1e-11));
    }
  }
}

TEST_CASE("YFS form factor: Kmax enters through the soft log only",
          "[YFS][FormFactor]") {
  YFS::YFS_Form_Factor &ff(FormFactor());
  for (const Pair &x : Pairs()) {
    INFO(x.name);
    const double p1p2(x.p1*x.p2);
    const double coef(alpha0/M_PI*(p1p2*ff.A(p1p2, x.m1, x.m2) - 1.));
    const double E(x.p1[0] + x.p2[0]);
    const double k1(1e-4*E), k2(0.03*E);
    CHECK(Real(x, k2, 1e-10) - Real(x, k1, 1e-10)
          == Approx(coef*std::log(k2*k2/(k1*k1))).epsilon(1e-12));
  }
}

TEST_CASE("YFS form factor: A is symmetric and agrees between its overloads",
          "[YFS][FormFactor]") {
  YFS::YFS_Form_Factor &ff(FormFactor());
  for (const Pair &x : Pairs()) {
    INFO(x.name);
    const double p1p2(x.p1*x.p2);
    CHECK(ff.A(p1p2, x.m1, x.m2) == Approx(ff.A(p1p2, x.m2, x.m1)).epsilon(1e-14));
    // the four-vector form takes the masses from p.Mass(): agree to its digits
    CHECK(ff.A(x.p1, x.p2) == Approx(ff.A(p1p2, x.m1, x.m2)).epsilon(1e-5));
  }
  // at threshold, p1.p2 = m1 m2, A is defined as 0
  CHECK(ff.A(0.105658*0.000510999, 0.105658, 0.000510999) == 0.);
}

TEST_CASE("YFS form factor: light back-to-back pair gives the standard YFS result",
          "[YFS][FormFactor]") {
  const double sqrts(91.1876), me(0.000510999);
  const Pair x(MakePair("e+e- at the Z pole", sqrts, me, me));
  const double L(std::log(sqrts*sqrts/(me*me)));
  const double gamma(2.*alpha0/M_PI*(L - 1.));
  for (double eps : {1e-4, 1e-3, 1e-2}) {
    INFO("eps = " << eps);
    const double kmax(0.5*eps*sqrts);
    const double yfs(gamma*std::log(eps) + 0.25*gamma
                     + alpha0/M_PI*(M_PI*M_PI/3. - 0.5));
    CHECK(Real(x, kmax, 1e-10) + Virtual(x, 1e-10) == Approx(yfs).margin(1e-8));
  }
}
