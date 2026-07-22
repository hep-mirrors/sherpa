#include <catch2/catch_all.hpp>

#include "AMISIC++/Perturbative/QCD_Processes.H"

using namespace AMISIC;

// Spin- and colour-averaged QCD 2->2 |M|^2 (with the strong coupling factored
// out, as the AMISIC processes store it). Pinned at one physical Mandelstam
// point s+t+u=0; the expected numbers are computed independently of the
// production formulas. qqbar->qqbar is additionally checked against the
// Ellis-Stirling-Webber closed form, which is what guards the interference
// sign (a wrong sign gives 7.13 instead of 8.10 here).
namespace {
  const double s =  100.;
  const double t =  -30.;
  const double u =  -70.;  // s+t+u = 0
}

TEST_CASE("AMISIC QCD 2->2 matrix elements reproduce the known averaged |M|^2",
          "[AMISIC::QCD_Processes]") {

  SECTION("qqbar->qqbar matches the ESW closed form (interference sign)") {
    // ESW: (4/9)[(s^2+u^2)/t^2 + (t^2+u^2)/s^2] - (8/27) u^2/(s t)
    const double esw = 4./9. * ((s*s+u*u)/(t*t) + (t*t+u*u)/(s*s))
                     - 8./27. * (u*u)/(s*t);
    qqbar_qqbar me;
    me.Calc(s, t, u);
    CHECK(me() == Catch::Approx(esw).epsilon(1e-9));
    CHECK(me() == Catch::Approx(8.0997531).epsilon(1e-6));
    // the pre-fix wrong sign would give 7.1318519 here
    CHECK(me() != Catch::Approx(7.1318519).epsilon(1e-6));
  }

  SECTION("distinct-flavour t-channel: q1q2->q1q2") {
    q1q2_q1q2 me;
    me.Calc(s, t, u);
    CHECK(me() == Catch::Approx(4./9. * (s*s+u*u)/(t*t)).epsilon(1e-9));
    CHECK(me() == Catch::Approx(7.3580247).epsilon(1e-6));
  }

  SECTION("distinct-flavour s-channel: q1q1bar->q2q2bar") {
    q1q1bar_q2q2bar me;
    me.Calc(s, t, u);
    CHECK(me() == Catch::Approx(4./9. * (t*t+u*u)/(s*s)).epsilon(1e-9));
    CHECK(me() == Catch::Approx(0.2577778).epsilon(1e-6));
  }

  SECTION("gluon and mixed channels are pinned at the same point") {
    gg_gg gg;       gg.Calc(s, t, u);
    CHECK(gg()       == Catch::Approx(25.155051).epsilon(1e-6));
    gg_qqbar ggq;   ggq.Calc(s, t, u);
    CHECK(ggq()      == Catch::Approx(0.2428175).epsilon(1e-6));
    qqbar_gg qgg;   qgg.Calc(s, t, u);
    CHECK(qgg()      == Catch::Approx(0.8633511).epsilon(1e-6));
    qg_qg qg;       qg.Calc(s, t, u);
    CHECK(qg()       == Catch::Approx(17.5015873).epsilon(1e-6));
    qq_qq qq;       qq.Calc(s, t, u);
    CHECK(qq()       == Catch::Approx(3.4678760).epsilon(1e-6));
  }

  SECTION("all matrix elements are positive at a physical point") {
    gg_gg gg;             gg.Calc(s,t,u);   CHECK(gg()  > 0.);
    gg_qqbar ggq;         ggq.Calc(s,t,u);  CHECK(ggq() > 0.);
    qqbar_gg qgg;         qgg.Calc(s,t,u);  CHECK(qgg() > 0.);
    qg_qg qg;             qg.Calc(s,t,u);   CHECK(qg()  > 0.);
    qq_qq qq;             qq.Calc(s,t,u);   CHECK(qq()  > 0.);
    qqbar_qqbar qqb;      qqb.Calc(s,t,u);  CHECK(qqb() > 0.);
    q1q2_q1q2 q12;        q12.Calc(s,t,u);  CHECK(q12() > 0.);
    q1q1bar_q2q2bar q11;  q11.Calc(s,t,u);  CHECK(q11() > 0.);
  }
}
